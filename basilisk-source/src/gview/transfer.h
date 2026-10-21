
#ifndef TRANSFER_H
#define TRANSFER_H

#include <stdbool.h>
#include <stddef.h>

/**
The type of message returned by transfer_receive().
*/
enum transfer_message_type {
  TRANSFER_MESSAGE_DATA,
  TRANSFER_MESSAGE_CLOSED
};

/**
A message received from the transfer.

For TRANSFER_MESSAGE_DATA, `data` is dynamically allocated and
belongs to the caller. The caller must call transfer_message_free()
when it no longer needs the message.

For TRANSFER_MESSAGE_CLOSED, `data` is NULL and `size` is zero.
*/
struct transfer_message {
  enum transfer_message_type type;
  void *data;
  size_t size;
  bool binary;
};

struct transfer;

/**
Open a connection to the specified host and port.

If `reconnect_ms` is zero, `transfer_open()` blocks until the
connection and WebSocket handshake have completed. It returns NULL
if the connection cannot be established.

If `reconnect_ms` is non-zero, `transfer_open()` does not wait for
the connection. The connection is established asynchronously by
`transfer_receive()`. If the connection is subsequently lost,
`transfer_receive()` will periodically attempt to reconnect, with
at least `reconnect_ms` milliseconds between attempts.

@param host Host name or IP address of the server.
@param port Port number of the server.
@param reconnect_ms Minimum interval between reconnection attempts,
or zero for a synchronous initial connection with no reconnection.

@return A new transfer object, or NULL if the synchronous connection
attempt fails.
*/
struct transfer * transfer_open (const char *host, int port,
                                 unsigned int reconnect_ms);

/**
Send a message to the server.

@param transfer Transfer object.
@param data Message data.
@param size Message size in bytes.
@param binary If true, send a binary WebSocket frame. Otherwise send
a text WebSocket frame.

@return 0 on success, -1 on error.
*/
int transfer_send (struct transfer *transfer,
                   const void *data, size_t size, bool binary);

/**
Poll for a received message.

This function is non-blocking. It returns NULL if no message is
currently available.

When reconnecting, this function also progresses the connection and
WebSocket handshake and performs reconnection attempts when the
reconnect timer expires.

@return A received message, or NULL if none is available.
*/
struct transfer_message * transfer_receive (struct transfer *transfer);

/**
Enable optional receive throttling.

The transfer layer monitors both the number of messages and the
total number of bytes waiting in its receive queue.

If either high-water limit is reached, a NOT_READY control message
is sent to the server.

A READY control message is sent when both low-water limits are
satisfied.

The control messages are an application-level convention and are
not part of the WebSocket protocol. Throttling is disabled by
default.

A high-water value of zero disables that particular limit. The
corresponding low-water value is then ignored.

@param transfer Transfer object.
@param high_bytes Queue size at which the server is asked to stop
sending, or zero to disable the byte limit.
@param low_bytes Queue size below which the byte limit is satisfied.
@param high_messages Number of queued messages at which the server
is asked to stop sending, or zero to disable the message limit.
@param low_messages Number of queued messages below which the
message limit is satisfied.

@return 0 on success, -1 on invalid arguments.
*/
int transfer_set_throttle (struct transfer *transfer,
                           size_t high_bytes, size_t low_bytes,
                           size_t high_messages, size_t low_messages);

/**
Disable receive throttling.

No further READY or NOT_READY control messages are generated.
*/
void transfer_disable_throttle (struct transfer *transfer);

/**
Free a message returned by transfer_receive().
*/
void transfer_message_free (struct transfer_message *message);

/**
Close the transfer and release all associated resources.
*/
int transfer_close (struct transfer *transfer);

#endif
