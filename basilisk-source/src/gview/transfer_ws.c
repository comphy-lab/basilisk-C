#define _POSIX_C_SOURCE 200809L

#include "transfer.h"

#include <arpa/inet.h>
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <time.h>
#include <unistd.h>
#include <netdb.h>

#define WS_RECV_BUFFER_SIZE 65536

#define TRANSFER_CONTROL_READY      1
#define TRANSFER_CONTROL_NOT_READY  0

static const unsigned char throttle_magic[] = {
  'T', 'R', 'A', 'N', 'S', 'F', 'E', 'R',
  '_', 'C', 'T', 'R', 'L'
};

enum connection_state {
  CONNECTION_DISCONNECTED,
  CONNECTION_CONNECTING,
  CONNECTION_HANDSHAKING,
  CONNECTION_CONNECTED
};

struct message_node {
  struct transfer_message *message;
  struct message_node *next;
};

struct transfer {
  char *host;
  int port;
  unsigned int reconnect_ms;

  int fd;
  enum connection_state state;

  struct timespec next_retry;

  char handshake_request[1024];
  size_t handshake_length;
  size_t handshake_sent;

  char handshake_response[4096];
  size_t handshake_received;

  /*
  Data received together with the HTTP handshake but belonging
  to the WebSocket stream.
  */
  unsigned char recv_buffer[WS_RECV_BUFFER_SIZE];
  size_t recv_start;
  size_t recv_end;

  pthread_t thread;
  bool thread_running;
  bool running;

  /*
  Protects the receive queue and connection state.
  */
  pthread_mutex_t mutex;

  /*
  Serialises WebSocket frame transmission. In particular, this
  prevents a throttling control frame from being interleaved with
  a frame sent by the application thread.
  */
  pthread_mutex_t send_mutex;

  struct message_node *queue_first;
  struct message_node *queue_last;

  size_t queue_bytes;
  size_t queue_messages;

  bool throttle;
  bool ready;

  size_t throttle_high_bytes;
  size_t throttle_low_bytes;

  size_t throttle_high_messages;
  size_t throttle_low_messages;
};


/*
 * ----------------------------------------------------------------------
 * Time and connection helpers
 * ----------------------------------------------------------------------
 */

static long long now_ms (void)
{
  struct timespec ts;

  clock_gettime (CLOCK_MONOTONIC, &ts);

  return (long long) ts.tv_sec * 1000 +
         ts.tv_nsec / 1000000;
}

static void set_next_retry (struct transfer *transfer)
{
  long long t;

  t = now_ms () + transfer->reconnect_ms;

  transfer->next_retry.tv_sec = t / 1000;
  transfer->next_retry.tv_nsec = (t % 1000) * 1000000;
}

static bool retry_due (struct transfer *transfer)
{
  long long retry;

  retry = (long long) transfer->next_retry.tv_sec * 1000 +
          transfer->next_retry.tv_nsec / 1000000;

  return now_ms () >= retry;
}

static int set_nonblocking (int fd)
{
  int flags;

  flags = fcntl (fd, F_GETFL, 0);

  if (flags < 0)
    return -1;

  return fcntl (fd, F_SETFL, flags | O_NONBLOCK);
}

static int set_blocking (int fd)
{
  int flags;

  flags = fcntl (fd, F_GETFL, 0);

  if (flags < 0)
    return -1;

  return fcntl (fd, F_SETFL, flags & ~O_NONBLOCK);
}

static void connection_failed (struct transfer *transfer)
{
  if (transfer->fd >= 0) {
    shutdown (transfer->fd, SHUT_RDWR);
    close (transfer->fd);
  }

  transfer->fd = -1;
  transfer->state = CONNECTION_DISCONNECTED;

  transfer->handshake_length = 0;
  transfer->handshake_sent = 0;
  transfer->handshake_received = 0;

  transfer->recv_start = 0;
  transfer->recv_end = 0;

  if (transfer->reconnect_ms > 0)
    set_next_retry (transfer);
}

static int start_connection (struct transfer *transfer)
{
  struct addrinfo hints, *result, *rp;
  char port[16];
  int fd = -1;
  int error;
  int connect_result;

  memset (&hints, 0, sizeof(hints));
  hints.ai_family = AF_INET;
  hints.ai_socktype = SOCK_STREAM;

  snprintf (port, sizeof(port), "%d", transfer->port);

  error = getaddrinfo (transfer->host, port, &hints, &result);
  if (error != 0)
    return -1;

  for (rp = result; rp; rp = rp->ai_next) {
    fd = socket (rp->ai_family, rp->ai_socktype, rp->ai_protocol);
    if (fd < 0)
      continue;

    if (set_nonblocking (fd) < 0) {
      close (fd);
      fd = -1;
      continue;
    }

    connect_result = connect (fd, rp->ai_addr, rp->ai_addrlen);

    if (connect_result == 0) {
      transfer->fd = fd;
      transfer->state = CONNECTION_HANDSHAKING;
      freeaddrinfo (result);
      return 0;
    }

    if (errno == EINPROGRESS) {
      transfer->fd = fd;
      transfer->state = CONNECTION_CONNECTING;
      freeaddrinfo (result);
      return 0;
    }

    close (fd);
    fd = -1;
  }

  freeaddrinfo (result);
  return -1;
}

static int finish_connect (struct transfer *transfer, int timeout)
{
  struct pollfd pfd;
  int error;
  socklen_t error_size;
  int result;

  pfd.fd = transfer->fd;
  pfd.events = POLLOUT;
  pfd.revents = 0;

  result = poll (&pfd, 1, timeout);

  if (result <= 0)
    return result;

  if (!(pfd.revents & (POLLOUT | POLLERR | POLLHUP)))
    return 0;

  error = 0;
  error_size = sizeof(error);

  if (getsockopt (transfer->fd, SOL_SOCKET, SO_ERROR,
                  &error, &error_size) < 0 || error != 0)
    return -1;

  transfer->state = CONNECTION_HANDSHAKING;

  return 1;
}


/*
 * ----------------------------------------------------------------------
 * WebSocket handshake
 * ----------------------------------------------------------------------
 */

static int start_handshake (struct transfer *transfer)
{
  int length;

  /*
   * The key is deliberately fixed. The client and server are under
   * our control and the server does not require validation of the
   * WebSocket key.
   */
  length = snprintf (transfer->handshake_request,
                     sizeof(transfer->handshake_request),
                     "GET / HTTP/1.1\r\n"
                     "Host: %s\r\n"
                     "Upgrade: websocket\r\n"
                     "Connection: Upgrade\r\n"
                     "Sec-WebSocket-Key: SGVsbG9XZWJTb2NrZXQ=\r\n"
                     "Sec-WebSocket-Version: 13\r\n"
                     "\r\n",
                     transfer->host);

  if (length < 0 ||
      (size_t) length >= sizeof(transfer->handshake_request))
    return -1;

  transfer->handshake_length = length;
  transfer->handshake_sent = 0;
  transfer->handshake_received = 0;

  transfer->recv_start = 0;
  transfer->recv_end = 0;

  return 0;
}

static int progress_handshake (struct transfer *transfer)
{
  ssize_t n;
  char buffer[1024];
  char *end;
  size_t header_size;
  size_t extra;

  while (transfer->handshake_sent < transfer->handshake_length) {
    n = send (transfer->fd,
              transfer->handshake_request +
              transfer->handshake_sent,
              transfer->handshake_length -
              transfer->handshake_sent,
              MSG_NOSIGNAL | MSG_DONTWAIT);

    if (n > 0) {
      transfer->handshake_sent += n;
      continue;
    }

    if (n < 0 && (errno == EAGAIN || errno == EWOULDBLOCK))
      return 0;

    return -1;
  }

  for (;;) {
    if (transfer->handshake_received >=
        sizeof(transfer->handshake_response) - 1)
      return -1;

    n = recv (transfer->fd, buffer, sizeof(buffer), MSG_DONTWAIT);

    if (n > 0) {
      if (transfer->handshake_received + (size_t) n >=
          sizeof(transfer->handshake_response))
        return -1;

      memcpy (transfer->handshake_response +
              transfer->handshake_received,
              buffer, n);

      transfer->handshake_received += n;

      transfer->handshake_response[
        transfer->handshake_received] = '\0';

      end = strstr (transfer->handshake_response, "\r\n\r\n");

      if (!end)
        continue;

      header_size = (size_t) (end + 4 -
                              transfer->handshake_response);

      if (strncmp (transfer->handshake_response,
                   "HTTP/1.1 101", 12) != 0)
        return -1;

      /*
       * Preserve any bytes following the HTTP headers. They already
       * belong to the WebSocket stream.
       */
      extra = transfer->handshake_received - header_size;

      if (extra > WS_RECV_BUFFER_SIZE)
        return -1;

      if (extra)
        memcpy (transfer->recv_buffer,
                transfer->handshake_response + header_size,
                extra);

      transfer->recv_start = 0;
      transfer->recv_end = extra;

      if (set_blocking (transfer->fd) < 0)
        return -1;

      transfer->state = CONNECTION_CONNECTED;

      return 1;
    }

    if (n == 0)
      return -1;

    if (errno == EAGAIN || errno == EWOULDBLOCK)
      return 0;

    return -1;
  }
}


/*
 * ----------------------------------------------------------------------
 * Receive queue
 * ----------------------------------------------------------------------
 */

static bool queue_push (struct transfer *transfer,
                        struct transfer_message *message);

static struct transfer_message * queue_pop (struct transfer *transfer);


/*
 * ----------------------------------------------------------------------
 * WebSocket transmission
 * ----------------------------------------------------------------------
 */

static int send_all (int fd, const void *buf, size_t n)
{
  const char *p = buf;

  while (n) {
    ssize_t r;

    r = send (fd, p, n, MSG_NOSIGNAL);

    if (r <= 0)
      return -1;

    p += r;
    n -= r;
  }

  return 0;
}

static int ws_send (struct transfer *transfer,
                    const void *data, size_t size, int opcode)
{
  uint8_t header[10];
  uint8_t mask[4];
  uint8_t buffer[4096];
  const uint8_t *src = data;
  size_t header_size;
  size_t n;
  size_t i;

  /*
   * A WebSocket client must mask frames sent to the server.
   *
   * This is deliberately not cryptographically random; the protocol
   * only requires a different masking key for each frame.
   */
  mask[0] = 0x12;
  mask[1] = 0x34;
  mask[2] = 0x56;
  mask[3] = 0x78;

  header[0] = 0x80 | opcode;

  if (size <= 125) {
    header[1] = 0x80 | size;
    header_size = 2;
  }
  else if (size <= 65535) {
    header[1] = 0x80 | 126;
    header[2] = size >> 8;
    header[3] = size;
    header_size = 4;
  }
  else {
    header[1] = 0x80 | 127;

    for (i = 0; i < 8; i++)
      header[2 + i] = size >> (56 - 8*i);

    header_size = 10;
  }

  if (send_all (transfer->fd, header, header_size) < 0 ||
      send_all (transfer->fd, mask, 4) < 0)
    return -1;

  while (size) {
    n = size < sizeof(buffer) ? size : sizeof(buffer);

    for (i = 0; i < n; i++)
      buffer[i] = src[i] ^ mask[i & 3];

    if (send_all (transfer->fd, buffer, n) < 0)
      return -1;

    src += n;
    size -= n;
  }

  return 0;
}

static int transfer_send_control (struct transfer *transfer,
                                  bool ready)
{
  unsigned char message[sizeof(throttle_magic) + 1];
  int result;

  memcpy (message, throttle_magic, sizeof(throttle_magic));

  message[sizeof(throttle_magic)] =
    ready ? TRANSFER_CONTROL_READY :
            TRANSFER_CONTROL_NOT_READY;

  pthread_mutex_lock (&transfer->send_mutex);

  pthread_mutex_lock (&transfer->mutex);

  if (transfer->state != CONNECTION_CONNECTED) {
    pthread_mutex_unlock (&transfer->mutex);
    pthread_mutex_unlock (&transfer->send_mutex);
    return -1;
  }

  pthread_mutex_unlock (&transfer->mutex);

  result = ws_send (transfer, message, sizeof(message), 0x2);

  pthread_mutex_unlock (&transfer->send_mutex);

  return result;
}


/*
 * ----------------------------------------------------------------------
 * Queue implementation
 * ----------------------------------------------------------------------
 */

static bool queue_push (struct transfer *transfer,
                        struct transfer_message *message)
{
  struct message_node *node;
  bool send_not_ready = false;

  node = malloc (sizeof(*node));

  if (!node) {
    transfer_message_free (message);
    return false;
  }

  node->message = message;
  node->next = NULL;

  pthread_mutex_lock (&transfer->mutex);

  if (transfer->queue_last)
    transfer->queue_last->next = node;
  else
    transfer->queue_first = node;

  transfer->queue_last = node;

  transfer->queue_bytes += message->size;
  transfer->queue_messages++;

  if (transfer->throttle &&
      transfer->ready &&
      ((transfer->throttle_high_bytes &&
        transfer->queue_bytes >= transfer->throttle_high_bytes) ||
       (transfer->throttle_high_messages &&
        transfer->queue_messages >=
        transfer->throttle_high_messages))) {
    transfer->ready = false;
    send_not_ready = true;
  }

  pthread_mutex_unlock (&transfer->mutex);

  if (send_not_ready)
    transfer_send_control (transfer, false);

  return true;
}

static struct transfer_message * queue_pop (struct transfer *transfer)
{
  struct message_node *node;
  struct transfer_message *message;
  bool send_ready = false;

  pthread_mutex_lock (&transfer->mutex);

  node = transfer->queue_first;

  if (node) {
    transfer->queue_first = node->next;

    if (!transfer->queue_first)
      transfer->queue_last = NULL;

    message = node->message;

    transfer->queue_bytes -= message->size;
    transfer->queue_messages--;

    if (transfer->throttle &&
        !transfer->ready &&
        (!transfer->throttle_low_bytes ||
         transfer->queue_bytes <= transfer->throttle_low_bytes) &&
        (!transfer->throttle_low_messages ||
         transfer->queue_messages <=
         transfer->throttle_low_messages)) {
      transfer->ready = true;
      send_ready = true;
    }
  }

  pthread_mutex_unlock (&transfer->mutex);

  if (!node)
    return NULL;

  free (node);

  if (send_ready)
    transfer_send_control (transfer, true);

  return message;
}


/*
 * ----------------------------------------------------------------------
 * WebSocket reception
 * ----------------------------------------------------------------------
 */

static ssize_t ws_buffer_recv (struct transfer *transfer,
                               void *buf, size_t size)
{
  size_t available;
  size_t n;

  available = transfer->recv_end - transfer->recv_start;

  if (available) {
    n = available < size ? available : size;

    memcpy (buf,
            transfer->recv_buffer + transfer->recv_start,
            n);

    transfer->recv_start += n;

    if (transfer->recv_start == transfer->recv_end) {
      transfer->recv_start = 0;
      transfer->recv_end = 0;
    }

    return n;
  }

  return recv (transfer->fd, buf, size, 0);
}

static int recv_all (struct transfer *transfer,
                     void *buf, size_t n)
{
  char *p = buf;

  while (n) {
    ssize_t r;

    r = ws_buffer_recv (transfer, p, n);

    if (r <= 0)
      return -1;

    p += r;
    n -= r;
  }

  return 0;
}

static int ws_receive (struct transfer *transfer)
{
  uint8_t header[2];
  uint8_t *payload = NULL;
  uint8_t opcode;
  uint64_t length;
  size_t payload_size;
  struct transfer_message *message;

  if (recv_all (transfer, header, 2) < 0)
    return -1;

  opcode = header[0] & 0x0f;
  length = header[1] & 0x7f;

  /*
   * The server sends unmasked frames.
   */
  if (header[1] & 0x80)
    return -1;

  if (length == 126) {
    uint8_t extended[2];

    if (recv_all (transfer, extended, 2) < 0)
      return -1;

    length = ((uint64_t) extended[0] << 8) |
             (uint64_t) extended[1];
  }
  else if (length == 127) {
    uint8_t extended[8];

    if (recv_all (transfer, extended, 8) < 0)
      return -1;

    length = ((uint64_t) extended[0] << 56) |
             ((uint64_t) extended[1] << 48) |
             ((uint64_t) extended[2] << 40) |
             ((uint64_t) extended[3] << 32) |
             ((uint64_t) extended[4] << 24) |
             ((uint64_t) extended[5] << 16) |
             ((uint64_t) extended[6] << 8) |
             (uint64_t) extended[7];
  }

  if (length > SIZE_MAX)
    return -1;

  payload_size = (size_t) length;

  if (payload_size) {
    payload = malloc (payload_size);

    if (!payload)
      return -1;

    if (recv_all (transfer, payload, payload_size) < 0) {
      free (payload);
      return -1;
    }
  }

  switch (opcode) {

  case 0x1:
  case 0x2:
    message = malloc (sizeof(*message));

    if (!message) {
      free (payload);
      return -1;
    }

    message->type = TRANSFER_MESSAGE_DATA;
    message->data = payload;
    message->size = payload_size;
    message->binary = opcode == 0x2;

    if (!queue_push (transfer, message))
      return -1;

    return 1;

  case 0x8:
    /*
     * The remote side closed the connection.
     */
    free (payload);
    return 0;

  case 0x9:
    /*
     * Ping. We do not expose it to the application.
     *
     * Respond with pong.
     */
    pthread_mutex_lock (&transfer->send_mutex);

    if (transfer->state == CONNECTION_CONNECTED)
      ws_send (transfer, payload, payload_size, 0xa);

    pthread_mutex_unlock (&transfer->send_mutex);

    free (payload);
    return 1;

  case 0xa:
    /*
     * Pong.
     */
    free (payload);
    return 1;

  default:
    free (payload);
    return 1;
  }
}


/*
 * ----------------------------------------------------------------------
 * Receiver thread
 * ----------------------------------------------------------------------
 */

static void queue_closed (struct transfer *transfer)
{
  struct transfer_message *message;

  message = malloc (sizeof(*message));

  if (!message)
    return;

  message->type = TRANSFER_MESSAGE_CLOSED;
  message->data = NULL;
  message->size = 0;
  message->binary = false;

  queue_push (transfer, message);
}

static void *receive_thread (void *arg)
{
  struct transfer *transfer = arg;
  int fd;
  int result;
  bool notify_closed = false;

  pthread_mutex_lock (&transfer->mutex);

  fd = transfer->fd;

  pthread_mutex_unlock (&transfer->mutex);

  for (;;) {
    result = ws_receive (transfer);

    if (result <= 0)
      break;
  }

  pthread_mutex_lock (&transfer->mutex);

  if (transfer->fd == fd) {
    close (transfer->fd);
    transfer->fd = -1;
  }

  transfer->state = CONNECTION_DISCONNECTED;
  transfer->thread_running = false;

  if (transfer->running) {
    notify_closed = true;

    if (transfer->reconnect_ms > 0)
      set_next_retry (transfer);
  }

  pthread_mutex_unlock (&transfer->mutex);

  if (notify_closed)
    queue_closed (transfer);

  return NULL;
}

static int start_receiver (struct transfer *transfer)
{
  if (pthread_create (&transfer->thread, NULL,
                      receive_thread, transfer) != 0)
    return -1;

  transfer->thread_running = true;

  return 0;
}


/*
 * ----------------------------------------------------------------------
 * Public interface
 * ----------------------------------------------------------------------
 */

struct transfer * transfer_open (const char *host, int port,
                                 unsigned int reconnect_ms)
{
  struct transfer *transfer;

  transfer = calloc (1, sizeof(*transfer));

  if (!transfer)
    return NULL;

  transfer->host = strdup (host);

  if (!transfer->host) {
    free (transfer);
    return NULL;
  }

  transfer->port = port;
  transfer->reconnect_ms = reconnect_ms;
  transfer->fd = -1;
  transfer->state = CONNECTION_DISCONNECTED;
  transfer->running = true;
  transfer->ready = true;

  if (pthread_mutex_init (&transfer->mutex, NULL) != 0) {
    free (transfer->host);
    free (transfer);
    return NULL;
  }

  if (pthread_mutex_init (&transfer->send_mutex, NULL) != 0) {
    pthread_mutex_destroy (&transfer->mutex);
    free (transfer->host);
    free (transfer);
    return NULL;
  }

  /*
   * With reconnect_ms == 0, establish the connection synchronously.
   */
  if (reconnect_ms == 0) {
    struct pollfd pfd;
    int result;

    if (start_connection (transfer) < 0)
      goto fail;

    if (transfer->state == CONNECTION_CONNECTING) {
      result = finish_connect (transfer, -1);

      if (result <= 0)
        goto fail;
    }

    if (start_handshake (transfer) < 0)
      goto fail;

    /*
     * progress_handshake() is non-blocking, so wait for the socket
     * whenever the handshake needs more data.
     */
    while (transfer->state == CONNECTION_HANDSHAKING) {
      pfd.fd = transfer->fd;
      pfd.events = POLLIN | POLLOUT;
      pfd.revents = 0;

      if (poll (&pfd, 1, -1) <= 0)
        goto fail;

      result = progress_handshake (transfer);

      if (result < 0)
        goto fail;
    }

    if (transfer->state != CONNECTION_CONNECTED)
      goto fail;

    /*
     * The asynchronous path starts the receiver immediately after
     * the handshake. Do the same here.
     */
    if (start_receiver (transfer) < 0)
      goto fail;
  }

  return transfer;

fail:
  if (transfer->fd >= 0) {
    shutdown (transfer->fd, SHUT_RDWR);
    close (transfer->fd);
  }

  pthread_mutex_destroy (&transfer->send_mutex);
  pthread_mutex_destroy (&transfer->mutex);

  free (transfer->host);
  free (transfer);

  return NULL;
}

struct transfer_message * transfer_receive (struct transfer *transfer)
{
  struct transfer_message *message;
  int result;

  /*
   * Messages already received take priority.
   */
  message = queue_pop (transfer);

  if (message)
    return message;

  pthread_mutex_lock (&transfer->mutex);

  switch (transfer->state) {

  case CONNECTION_DISCONNECTED:

    /*
     * With reconnect_ms == 0, transfer_open() made the only
     * connection attempt. Never reconnect here.
     */
    if (transfer->reconnect_ms == 0) {
      pthread_mutex_unlock (&transfer->mutex);
      return NULL;
    }

    if (!retry_due (transfer)) {
      pthread_mutex_unlock (&transfer->mutex);
      return NULL;
    }

    set_next_retry (transfer);

    result = start_connection (transfer);

    if (result < 0) {
      pthread_mutex_unlock (&transfer->mutex);
      return NULL;
    }

    if (transfer->state == CONNECTION_HANDSHAKING) {
      if (start_handshake (transfer) < 0) {
        connection_failed (transfer);
        pthread_mutex_unlock (&transfer->mutex);
        return NULL;
      }
    }

    break;

  case CONNECTION_CONNECTING:

    result = finish_connect (transfer, 0);

    if (result < 0) {
      connection_failed (transfer);
      pthread_mutex_unlock (&transfer->mutex);
      return NULL;
    }

    if (transfer->state == CONNECTION_HANDSHAKING) {
      if (start_handshake (transfer) < 0) {
        connection_failed (transfer);
        pthread_mutex_unlock (&transfer->mutex);
        return NULL;
      }
    }

    break;

  case CONNECTION_HANDSHAKING:

    result = progress_handshake (transfer);

    if (result < 0) {
      connection_failed (transfer);
      pthread_mutex_unlock (&transfer->mutex);
      return NULL;
    }

    if (result > 0) {
      if (start_receiver (transfer) < 0) {
        connection_failed (transfer);
        pthread_mutex_unlock (&transfer->mutex);
        return NULL;
      }
    }

    break;

  case CONNECTION_CONNECTED:
    break;
  }

  pthread_mutex_unlock (&transfer->mutex);

  return queue_pop (transfer);
}

int transfer_send (struct transfer *transfer,
                   const void *data, size_t size, bool binary)
{
  int result;

  pthread_mutex_lock (&transfer->send_mutex);

  pthread_mutex_lock (&transfer->mutex);

  if (transfer->state != CONNECTION_CONNECTED) {
    pthread_mutex_unlock (&transfer->mutex);
    pthread_mutex_unlock (&transfer->send_mutex);
    return -1;
  }

  pthread_mutex_unlock (&transfer->mutex);

  result = ws_send (transfer, data, size,
                    binary ? 0x2 : 0x1);

  pthread_mutex_unlock (&transfer->send_mutex);

  if (result < 0) {
    pthread_mutex_lock (&transfer->mutex);
    connection_failed (transfer);
    pthread_mutex_unlock (&transfer->mutex);
  }

  return result;
}

int transfer_set_throttle (struct transfer *transfer,
                           size_t high_bytes, size_t low_bytes,
                           size_t high_messages, size_t low_messages)
{
  bool send_not_ready = false;

  if (!transfer)
    return -1;

  if (high_bytes && low_bytes >= high_bytes)
    return -1;

  if (high_messages && low_messages >= high_messages)
    return -1;

  pthread_mutex_lock (&transfer->mutex);

  transfer->throttle_high_bytes = high_bytes;
  transfer->throttle_low_bytes = low_bytes;

  transfer->throttle_high_messages = high_messages;
  transfer->throttle_low_messages = low_messages;

  transfer->throttle = true;
  transfer->ready = true;

  if ((high_bytes &&
       transfer->queue_bytes >= high_bytes) ||
      (high_messages &&
       transfer->queue_messages >= high_messages)) {
    transfer->ready = false;
    send_not_ready = true;
  }

  pthread_mutex_unlock (&transfer->mutex);

  if (send_not_ready)
    return transfer_send_control (transfer, false);

  return 0;
}

void transfer_disable_throttle (struct transfer *transfer)
{
  if (!transfer)
    return;

  pthread_mutex_lock (&transfer->mutex);

  transfer->throttle = false;
  transfer->ready = true;

  pthread_mutex_unlock (&transfer->mutex);
}

void transfer_message_free (struct transfer_message *message)
{
  if (!message)
    return;

  free (message->data);
  free (message);
}

int transfer_close (struct transfer *transfer)
{
  struct transfer_message *message;
  pthread_t thread;
  bool thread_running;

  if (!transfer)
    return 0;

  pthread_mutex_lock (&transfer->mutex);

  transfer->running = false;

  thread_running = transfer->thread_running;
  thread = transfer->thread;

  if (transfer->fd >= 0)
    shutdown (transfer->fd, SHUT_RDWR);

  pthread_mutex_unlock (&transfer->mutex);

  if (thread_running)
    pthread_join (thread, NULL);

  pthread_mutex_lock (&transfer->mutex);

  if (transfer->fd >= 0) {
    close (transfer->fd);
    transfer->fd = -1;
  }

  pthread_mutex_unlock (&transfer->mutex);

  while ((message = queue_pop (transfer)))
    transfer_message_free (message);

  pthread_mutex_destroy (&transfer->send_mutex);
  pthread_mutex_destroy (&transfer->mutex);

  free (transfer->host);
  free (transfer);

  return 0;
}
