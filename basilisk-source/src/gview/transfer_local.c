/**
# Implementation of "pseudo-websocket" local transfers

The server and client run on separate threads, but in the same process
i.e. the communication is purely "local" to the process. */

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>

/**
This implements both the client interface in [transfer.h]() used by
[display.c]() and the server interface in [../wsServer/include/ws.h]()
used by [../display.h](). */

#include "transfer.h"
#include "../wsServer/include/ws.h"

#define LOCAL_FD 1

struct local_message {
  struct local_message *next;
  char *data;
  size_t size;
  int type;
};

struct local_queue {
  struct local_message *first;
  struct local_message *last;
  pthread_mutex_t mutex;
};

struct local_transfer {
  bool closed;
};

static struct {
  bool server_open;
  bool closed;

  struct local_queue to_gview;
  struct local_queue to_server;

  /**
  Temporary buffer used by ws_sendframe_init()/ws_send(). */
  
  char *send_data;
  size_t send_size;
  size_t send_pos;
  int send_type;
  pthread_mutex_t send_mutex;
} local = {
  .server_open = false,
  .closed = false,

  .to_gview = {
    .first = NULL,
    .last = NULL,
    .mutex = PTHREAD_MUTEX_INITIALIZER
  },

  .to_server = {
    .first = NULL,
    .last = NULL,
    .mutex = PTHREAD_MUTEX_INITIALIZER
  },

  .send_data = NULL,
  .send_size = 0,
  .send_pos = 0,
  .send_type = 0,
  .send_mutex = PTHREAD_MUTEX_INITIALIZER
};

/**
Queue a copy of a message.

Ownership of the queued copy belongs to the queue until it is returned
to the receiving side. */

static int queue_put (struct local_queue *queue,
                      const void *data, size_t size, int type)
{
  struct local_message *message;

  message = malloc (sizeof(*message));
  if (!message)
    return -1;

  message->data = NULL;
  if (size) {
    message->data = malloc (size);
    if (!message->data) {
      free (message);
      return -1;
    }
    memcpy (message->data, data, size);
  }

  message->size = size;
  message->type = type;
  message->next = NULL;

  pthread_mutex_lock (&queue->mutex);

  if (queue->last)
    queue->last->next = message;
  else
    queue->first = message;

  queue->last = message;

  pthread_mutex_unlock (&queue->mutex);

  return 0;
}

static struct local_message * queue_get (struct local_queue *queue)
{
  struct local_message *message;

  pthread_mutex_lock (&queue->mutex);

  message = queue->first;
  if (message) {
    queue->first = message->next;
    if (!queue->first)
      queue->last = NULL;
  }

  pthread_mutex_unlock (&queue->mutex);

  return message;
}

static void queue_clear (struct local_queue *queue)
{
  struct local_message *message;

  while ((message = queue_get (queue))) {
    free (message->data);
    free (message);
  }
}

/**
# wsServer interface */

int ws_socket_open (uint16_t port)
{
  (void) port;

  if (local.server_open)
    return LOCAL_FD;

  local.server_open = true;
  local.closed = false;

  return LOCAL_FD;
}

int ws_sendframe_init (int fd, ssize_t size,
                       bool broadcast, int type)
{
  (void) broadcast;

  if (fd != LOCAL_FD || size < 0 || local.closed)
    return -1;

  pthread_mutex_lock (&local.send_mutex);

  if (local.send_data) {
    pthread_mutex_unlock (&local.send_mutex);
    return -1;
  }

  local.send_data = NULL;

  if (size > 0) {
    local.send_data = malloc ((size_t) size);
    if (!local.send_data) {
      pthread_mutex_unlock (&local.send_mutex);
      return -1;
    }
  }

  local.send_size = (size_t) size;
  local.send_pos = 0;
  local.send_type = type;

  /**
  An empty frame is complete immediately. */
  
  if (size == 0) {
    if (queue_put (&local.to_gview, NULL, 0, type) < 0) {
      free (local.send_data);
      local.send_data = NULL;
      local.send_size = 0;
      local.send_type = 0;
      pthread_mutex_unlock (&local.send_mutex);
      return -1;
    }

    local.send_data = NULL;
    local.send_size = 0;
    local.send_pos = 0;
    local.send_type = 0;
  }

  pthread_mutex_unlock (&local.send_mutex);

  return 0;
}

ssize_t ws_send (int fd, const void *buf, size_t len)
{
  int status = 0;

  if (fd != LOCAL_FD || local.closed)
    return -1;

  pthread_mutex_lock (&local.send_mutex);

  if (!local.send_data && local.send_size != 0) {
    pthread_mutex_unlock (&local.send_mutex);
    return -1;
  }

  if (local.send_pos + len > local.send_size) {
    pthread_mutex_unlock (&local.send_mutex);
    return -1;
  }

  if (len)
    memcpy (local.send_data + local.send_pos, buf, len);

  local.send_pos += len;

  /**
  ws_sendframe_init() specifies the complete frame size.
  Publish the message once all of it has been supplied. */
  
  if (local.send_pos == local.send_size) {
    if (queue_put (&local.to_gview,
                   local.send_data,
                   local.send_size,
                   local.send_type) < 0)
      status = -1;

    free (local.send_data);
    local.send_data = NULL;
    local.send_size = 0;
    local.send_pos = 0;
    local.send_type = 0;
  }

  pthread_mutex_unlock (&local.send_mutex);

  return status < 0 ? -1 : (ssize_t) len;
}

int ws_sendframe (int fd, const char *msg, ssize_t size,
                  bool broadcast, int type)
{
  (void) broadcast;

  if (fd != LOCAL_FD || size < 0 || local.closed)
    return -1;

  return queue_put (&local.to_gview, msg, (size_t) size, type);
}

int ws_sendframe_txt (int fd, const char *msg, bool broadcast)
{
  return ws_sendframe (fd, msg, strlen (msg),
                       broadcast, WS_FR_OP_TXT);
}

struct ws_message * ws_socket_poll (int fd, int timeout)
{
  struct local_message *message;
  struct ws_message *messages;

  (void) timeout;

  if (fd != LOCAL_FD)
    return NULL;

  message = queue_get (&local.to_server);
  if (!message)
    return NULL;

  /**
  ws_socket_poll() returns a NULL-terminated array in the
  wsServer convention: the final entry has fd < 0. */
  
  messages = malloc (2*sizeof(*messages));
  if (!messages) {
    free (message->data);
    free (message);
    return NULL;
  }

  messages[0].fd = LOCAL_FD;
  messages[0].msg = message->data;
  messages[0].size = message->size;
  messages[0].type = message->type;

  messages[1].fd = -1;
  messages[1].msg = NULL;
  messages[1].size = 0;
  messages[1].type = 0;

  free (message);

  return messages;
}

/**
# transfer interface */

struct transfer * transfer_open (const char *host, int port,
                                 unsigned int reconnect_ms)
{
  struct local_transfer *transfer;

  (void) host;
  (void) port;
  (void) reconnect_ms;

  /**
  In the local implementation the server must already have been opened
  by ws_socket_open(). There is no asynchronous connection or
  reconnection. */
  
  if (!local.server_open || local.closed)
    return NULL;

  transfer = calloc (1, sizeof(*transfer));
  if (!transfer)
    return NULL;

  if (queue_put (&local.to_server, NULL, 0, WS_FR_OP_OPEN) < 0) {
    free (transfer);
    return NULL;
  }
  
  return (struct transfer *) transfer;
}

int transfer_send (struct transfer *transfer,
                   const void *data, size_t size, bool binary)
{
  struct local_transfer *t = (struct local_transfer *) transfer;
  int type;

  if (!t || t->closed || local.closed)
    return -1;

  type = binary ? WS_FR_OP_BIN : WS_FR_OP_TXT;

  return queue_put (&local.to_server, data, size, type);
}

struct transfer_message * transfer_receive (struct transfer *transfer)
{
  struct local_transfer *t = (struct local_transfer *) transfer;
  struct local_message *message;
  struct transfer_message *result;

  if (!t || t->closed)
    return NULL;

  message = queue_get (&local.to_gview);
  if (!message)
    return NULL;

  result = malloc (sizeof(*result));
  if (!result) {
    free (message->data);
    free (message);
    return NULL;
  }

  if (message->type == WS_FR_OP_CLSE) {
    result->type = TRANSFER_MESSAGE_CLOSED;
    result->data = NULL;
    result->size = 0;
    result->binary = false;

    free (message->data);
  }
  else {
    result->type = TRANSFER_MESSAGE_DATA;
    result->data = message->data;
    result->size = message->size;
    result->binary = message->type == WS_FR_OP_BIN;
  }

  free (message);

  return result;
}

int transfer_set_throttle (struct transfer *transfer,
                           size_t high_bytes, size_t low_bytes,
                           size_t high_messages, size_t low_messages)
{
  (void) transfer;
  (void) high_bytes;
  (void) low_bytes;
  (void) high_messages;
  (void) low_messages;

  /**
  There is no network peer to throttle in the local implementation. */
  
  return 0;
}

void transfer_disable_throttle (struct transfer *transfer)
{
  (void) transfer;
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
  struct local_transfer *t = (struct local_transfer *) transfer;

  if (!t)
    return -1;

  /**
  There is only one connection and its lifetime ends with the
  program, so no reconnection or intermediate connection state
  is required. */
  
  t->closed = true;
  local.closed = true;

  pthread_mutex_lock (&local.send_mutex);

  free (local.send_data);
  local.send_data = NULL;
  local.send_size = 0;
  local.send_pos = 0;
  local.send_type = 0;

  pthread_mutex_unlock (&local.send_mutex);

  queue_clear (&local.to_gview);
  queue_clear (&local.to_server);

  free (t);

  return 0;
}
