#include "transfer.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>

int main (void)
{
  struct transfer *transfer;
  struct transfer_message *message;
  char line[4096];
  size_t size;

  transfer = transfer_open ("loclahost", 7100, 0);
  usleep (1000);
 
  if (!transfer) {
    fprintf (stderr, "Could not connect\n");
    return 1;
  }

  bool transfer_ready = true;
  for (;;) {
    
    /**
    Poll for received messages. */
    
    while ((message = transfer_receive (transfer))) {
      if (message->type == TRANSFER_MESSAGE_CLOSED) {
        printf ("Connection closed\n");
        transfer_message_free (message);
        transfer_close (transfer);
        return 0;
      }

      if (message->binary)
        printf ("Received binary message: %zu bytes\n",
                message->size);
      else
        printf ("Received: %.*s\n",
                (int) message->size,
                (char *) message->data);

      transfer_message_free (message);
      transfer_ready = true;
    }

    /**
    Do whatever other work the application needs to do. */
    
    usleep (1000);

    /**
    Example: read a request from stdin and send it.
    
    This is just illustrative; in your actual application this
    would presumably be part of the normal polling/work loop. */
    
    if (transfer_ready && fgets (line, sizeof(line), stdin)) {
      size = strlen (line);

      if (size && line[size - 1] == '\n')
        size--;

      if (transfer_send (transfer, line, size, false) < 0)
        break;
      transfer_ready = false;
    }
  }

  transfer_close (transfer);

  return 0;
}
