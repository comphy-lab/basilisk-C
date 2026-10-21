/**
# Emscripten WebSocket implementation for Basilisk gview */

#include "transfer.h"

#include <emscripten.h>
#include <emscripten/html5.h>
#include <emscripten/websocket.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TRANSFER_CONTROL_READY      1
#define TRANSFER_CONTROL_NOT_READY  0

static const unsigned char throttle_magic[] = {
  'T', 'R', 'A', 'N', 'S', 'F', 'E', 'R',
  '_', 'C', 'T', 'R', 'L'
};

enum connection_state {
  CONNECTION_DISCONNECTED,
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

  EMSCRIPTEN_WEBSOCKET_T socket;
  enum connection_state state;

  long reconnect_timer;
  bool closing;

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

static bool queue_push (struct transfer *t, struct transfer_message *msg);
static struct transfer_message *queue_pop (struct transfer *t);

static void reconnect_timer_fired (void *userData);

static void schedule_reconnect (struct transfer *t)
{
  if (t->closing ||
      !t->reconnect_ms ||
      t->state == CONNECTION_CONNECTED ||
      t->reconnect_timer)
    return;

  t->reconnect_timer =
    emscripten_set_timeout (reconnect_timer_fired,
                            t->reconnect_ms,
                            t);
}

static void cancel_reconnect (struct transfer *t)
{
  if (t->reconnect_timer) {
    emscripten_clear_timeout (t->reconnect_timer);
    t->reconnect_timer = 0;
  }
}

static bool connect (struct transfer *t);

static void reconnect_timer_fired (void *userData)
{
  struct transfer *t = userData;

  /**
  The timeout has fired, so it is no longer pending. */
  
  t->reconnect_timer = 0;

  if (t->closing || t->state == CONNECTION_CONNECTED)
    return;

  if (!connect (t))
    schedule_reconnect (t);
}

static int send_control (struct transfer *t, bool ready)
{
  unsigned char msg[sizeof (throttle_magic) + 1];

  memcpy (msg, throttle_magic, sizeof (throttle_magic));

  msg[sizeof (throttle_magic)] =
    ready ? TRANSFER_CONTROL_READY : TRANSFER_CONTROL_NOT_READY;

  if (t->state != CONNECTION_CONNECTED)
    return -1;

  return emscripten_websocket_send_binary
    (t->socket, msg, sizeof (msg)) == EMSCRIPTEN_RESULT_SUCCESS ? 0 : -1;
}

static EM_BOOL on_open (int eventType,
                        const EmscriptenWebSocketOpenEvent *e,
                        void *userData)
{
  struct transfer *t = userData;

  if (t->closing)
    return EM_TRUE;

  t->state = CONNECTION_CONNECTED;

  /**
  There should normally not be a timer at this point, but cancel it
  defensively in case the connection was established before a
  previously scheduled timer fired. */
  
  cancel_reconnect (t);

  /**
  If throttling is enabled and we are currently not ready, make sure
  the server receives the current state after reconnecting. */
  
  if (t->throttle && !t->ready)
    send_control (t, false);

  return EM_TRUE;
}


static void report_disconnect (struct transfer *t)
{
  /**
  Avoid generating multiple CLOSED messages for the same
  connection. In particular, a WebSocket error can be followed by a
  close event. */
  
  if (t->state == CONNECTION_DISCONNECTED)
    return;

  t->state = CONNECTION_DISCONNECTED;

  struct transfer_message *m = malloc (sizeof (*m));
  if (!m)
    return;

  m->type = TRANSFER_MESSAGE_CLOSED;
  m->data = NULL;
  m->size = 0;
  m->binary = false;

  queue_push (t, m);
}


static EM_BOOL on_error (int eventType,
                         const EmscriptenWebSocketErrorEvent *e,
                         void *userData)
{
  struct transfer *t = userData;

  if (t->closing)
    return EM_TRUE;

  /**
  Do not reconnect here. The browser will normally follow the error
  event with a close event, and reconnecting from on_close avoids
  creating a second WebSocket while the old one still exists.*/
  
  report_disconnect (t);

  return EM_TRUE;
}


static EM_BOOL on_close (int eventType,
                         const EmscriptenWebSocketCloseEvent *e,
                         void *userData)
{
  struct transfer *t = userData;

  if (t->closing)
    return EM_TRUE;

  /**
  The socket handle belongs to the WebSocket which just closed. */
  
  t->socket = 0;

  report_disconnect (t);

  /**
  Start a new connection after the requested delay. */
  
  schedule_reconnect (t);

  return EM_TRUE;
}


static EM_BOOL on_message (int eventType,
                           const EmscriptenWebSocketMessageEvent *e,
                           void *userData)
{
  struct transfer *t = userData;

  if (t->closing)
    return EM_TRUE;

  struct transfer_message *m = malloc (sizeof (*m));
  if (!m)
    return EM_TRUE;

  m->type = TRANSFER_MESSAGE_DATA;
  m->size = e->numBytes;
  m->binary = !e->isText;

  /**
  malloc(0) is allowed to return NULL, so don't treat NULL as
  an allocation failure for a zero-length message.*/
  
  m->data = e->numBytes ? malloc (e->numBytes) : NULL;
  if (e->numBytes) {
    if (!m->data) {
      free (m);
      return EM_TRUE;
    }
    memcpy (m->data, e->data, e->numBytes);
  }
    
  queue_push (t, m);

  return EM_TRUE;
}

static bool queue_push (struct transfer *t,
                        struct transfer_message *msg)
{
  struct message_node *node = malloc (sizeof (*node));

  if (!node) {
    transfer_message_free (msg);
    return false;
  }

  node->message = msg;
  node->next = NULL;

  if (t->queue_last)
    t->queue_last->next = node;
  else
    t->queue_first = node;

  t->queue_last = node;

  t->queue_bytes += msg->size;
  t->queue_messages++;

  if (t->throttle && t->ready &&
      ((t->throttle_high_bytes &&
        t->queue_bytes >= t->throttle_high_bytes) ||
       (t->throttle_high_messages &&
        t->queue_messages >= t->throttle_high_messages))) {

    t->ready = false;
    send_control (t, false);
  }

  return true;
}


static struct transfer_message *queue_pop (struct transfer *t)
{
  struct message_node *node = t->queue_first;

  if (!node)
    return NULL;

  t->queue_first = node->next;

  if (!t->queue_first)
    t->queue_last = NULL;

  struct transfer_message *msg = node->message;

  t->queue_bytes -= msg->size;
  t->queue_messages--;

  if (t->throttle && !t->ready &&
      (!t->throttle_low_bytes ||
       t->queue_bytes <= t->throttle_low_bytes) &&
      (!t->throttle_low_messages ||
       t->queue_messages <= t->throttle_low_messages)) {

    t->ready = true;
    send_control (t, true);
  }

  free (node);

  return msg;
}


static bool connect (struct transfer *t)
{
  if (t->closing)
    return false;

  /**
  Don't create another socket if one is already active. */
  
  if (t->socket)
    return false;

  char url[256];

  snprintf (url, sizeof (url), "ws://%s:%d", t->host, t->port);

  EmscriptenWebSocketCreateAttributes attr;

  emscripten_websocket_init_create_attributes (&attr);

  attr.url = url;

  t->socket = emscripten_websocket_new (&attr);

  if (!t->socket)
    return false;

  /**
  Register all callbacks before returning. */
  
  if (emscripten_websocket_set_onopen_callback
      (t->socket, t, on_open) != EMSCRIPTEN_RESULT_SUCCESS ||
      emscripten_websocket_set_onerror_callback
      (t->socket, t, on_error) != EMSCRIPTEN_RESULT_SUCCESS ||
      emscripten_websocket_set_onclose_callback
      (t->socket, t, on_close) != EMSCRIPTEN_RESULT_SUCCESS ||
      emscripten_websocket_set_onmessage_callback
      (t->socket, t, on_message) != EMSCRIPTEN_RESULT_SUCCESS) {

    emscripten_websocket_close (t->socket, 1000, NULL);
    emscripten_websocket_delete (t->socket);
    t->socket = 0;

    return false;
  }

  return true;
}


struct transfer *transfer_open (const char *host,
                                int port,
                                unsigned int reconnect_ms)
{
  struct transfer *t = calloc (1, sizeof (*t));

  if (!t)
    return NULL;

  t->host = strdup (host);

  if (!t->host) {
    free (t);
    return NULL;
  }

  t->port = port;
  t->reconnect_ms = reconnect_ms;

  t->socket = 0;
  t->state = CONNECTION_DISCONNECTED;
  t->reconnect_timer = 0;
  t->closing = false;

  t->ready = true;

  /**
  In a browser, creating the WebSocket is asynchronous.
  
  emscripten_websocket_new() only tells us whether the WebSocket
  object could be created; it does not mean that the server has
  accepted the connection yet.*/
  
  if (!connect (t)) {
    if (!reconnect_ms) {
      free (t->host);
      free (t);
      return NULL;
    }

    schedule_reconnect (t);
  }

  return t;
}


struct transfer_message *transfer_receive (struct transfer *t)
{
  if (!t)
    return NULL;

  return queue_pop (t);
}


int transfer_send (struct transfer *t,
                   const void *data,
                   size_t size,
                   bool binary)
{
  if (!t || t->closing ||
      t->state != CONNECTION_CONNECTED)
    return -1;

  if (binary) {
    return emscripten_websocket_send_binary
      (t->socket, (void *) data, size) ==
      EMSCRIPTEN_RESULT_SUCCESS ? 0 : -1;
  }

  /**
  The text API takes a NUL-terminated string, so size is not
  representable here. This is the same limitation as the native
  WebSocket text interface if it takes a C string. */
  
  (void) size;

  return emscripten_websocket_send_utf8_text
    (t->socket, data) ==
    EMSCRIPTEN_RESULT_SUCCESS ? 0 : -1;
}


int transfer_set_throttle (struct transfer *t,
                           size_t high_bytes,
                           size_t low_bytes,
                           size_t high_messages,
                           size_t low_messages)
{
  if (!t)
    return -1;

  if (high_bytes && low_bytes >= high_bytes)
    return -1;

  if (high_messages && low_messages >= high_messages)
    return -1;

  t->throttle_high_bytes = high_bytes;
  t->throttle_low_bytes = low_bytes;

  t->throttle_high_messages = high_messages;
  t->throttle_low_messages = low_messages;

  t->throttle = true;
  t->ready = true;

  if ((high_bytes && t->queue_bytes >= high_bytes) ||
      (high_messages && t->queue_messages >= high_messages)) {

    t->ready = false;

    /**
    If we are disconnected, this simply has no effect. The state will
    be sent by on_open() after reconnect. */
    
    send_control (t, false);
  }

  return 0;
}


void transfer_disable_throttle (struct transfer *t)
{
  if (!t)
    return;

  t->throttle = false;
  t->ready = true;
}


void transfer_message_free (struct transfer_message *msg)
{
  if (!msg)
    return;

  free (msg->data);
  free (msg);
}

int transfer_close (struct transfer *t)
{
  if (!t)
    return 0;

  t->closing = true;

  cancel_reconnect (t);

  if (t->socket) {
    EMSCRIPTEN_WEBSOCKET_T socket = t->socket;

    t->socket = 0;

    emscripten_websocket_close (socket, 1000, NULL);
    emscripten_websocket_delete (socket);
  }

  struct message_node *node = t->queue_first;

  while (node) {
    struct message_node *next = node->next;

    transfer_message_free (node->message);
    free (node);

    node = next;
  }

  free (t->host);
  free (t);

  return 0;
}
