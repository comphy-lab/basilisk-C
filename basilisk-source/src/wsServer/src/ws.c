/*
 * Copyright (C) 2016-2020  Davidson Francis <davidsondfgl@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#define _POSIX_C_SOURCE 200809L
#include <arpa/inet.h>
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <unistd.h>

#include <ws.h>

/**
 * @dir src/
 * @brief wsServer source code
 *
 * @file ws.c
 * @brief wsServer main routines.
 */

/**
 * @brief WebSocket frame data
 */
struct ws_frame_data
{
	/**
	 * @brief Frame read.
	 */
	unsigned char frm[MESSAGE_LENGTH];
	/**
	 * @brief Processed message at the moment.
	 */
	unsigned char *msg;
	/**
	 * @brief Control frame payload
	 */
	unsigned char msg_ctrl[125];
	/**
	 * @brief Current byte position.
	 */
	size_t cur_pos;
	/**
	 * @brief Amount of read bytes.
	 */
	size_t amt_read;
	/**
	 * @brief Frame type, like text or binary.
	 */
	int frame_type;
	/**
	 * @brief Frame size.
	 */
	size_t frame_size;
	/**
	 * @brief Error flag, set when a read was not possible.
	 */
	int error;
	/**
	 * @brief Client socket file descriptor.
	 */
	int sock;
};

/**
 * @brief Client socks.
 */
struct ws_connection
{
	int client_sock; /**< Client socket FD.        */
	int state;       /**< WebSocket current state. */
	struct ws_frame_data frame_data;
};

/**
 * @brief Clients list.
 */
struct ws_connection client_socks[MAX_CLIENTS];

/**
 * @brief Issues an error message and aborts the program.
 *
 * @param s Error message.
 */
#define panic(s)   \
	do             \
	{              \
		perror(s); \
		exit(-1);  \
	} while (0);

/**
 * @brief For a given client @p fd, returns its
 * client index if exists, or -1 otherwise.
 *
 * @param fd Client fd.
 *
 * @return Return the client index or -1 if invalid
 * fd.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int get_client_index(int fd)
{
	int i;
	for (i = 0; i < MAX_CLIENTS; i++)
		if (client_socks[i].client_sock == fd)
			break;
	return (i == MAX_CLIENTS ? -1 : i);
}

/**
 * @brief Returns the current client state for a given
 * client @p idx.
 *
 * @param idx Client index.
 *
 * @return Returns the client state, -1 otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int get_client_state(int idx)
{
	int state;

	if (idx < 0 || idx >= MAX_CLIENTS)
		return (-1);

	state = client_socks[idx].state;
	return (state);
}

/**
 * @brief Set a state @p state to the client index
 * @p idx.
 *
 * @param idx Client index.
 * @param state State to be set.
 *
 * @return Returns 0 if success, -1 otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int set_client_state(int idx, int state)
{
	if (idx < 0 || idx >= MAX_CLIENTS)
		return (-1);

	if (state < 0 || state > 3)
		return (-1);

	client_socks[idx].state = state;
	return (0);
}

#if 0 // fixme
/**
 * @brief Close time-out thread.
 *
 * For a given client, this routine sleeps until
 * TIMEOUT_MS and closes the connection or returns
 * sooner if already closed connection.
 *
 * @param p ws_connection Structure Pointer.
 *
 * @return Always NULL.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static void *close_timeout(void *p)
{
	struct ws_connection *conn = p;
	struct timespec ts;

	pthread_mutex_lock(&conn->mtx_state);

	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_nsec += MS_TO_NS(TIMEOUT_MS);

	/* Normalize the time. */
	while (ts.tv_nsec >= 1000000000)
	{
		ts.tv_sec++;
		ts.tv_nsec -= 1000000000;
	}

	while (conn->state != WS_STATE_CLOSED &&
		   pthread_cond_timedwait(&conn->cnd_state_close, &conn->mtx_state, &ts) !=
			   ETIMEDOUT)
		;

	/* If already closed. */
	if (conn->state == WS_STATE_CLOSED)
		goto quit;

	DEBUG("Timer expired, closing client %d\n", conn->client_sock);

	shutdown(conn->client_sock, SHUT_RDWR);
	close(conn->client_sock);
	conn->client_sock = -1;
	conn->state       = WS_STATE_CLOSED;
quit:
	pthread_mutex_unlock(&conn->mtx_state);
	return (NULL);
}

/**
 * @brief For a valid client index @p idx, starts
 * the timeout thread and set the current state
 * to 'CLOSING'.
 *
 * @param idx Client index.
 *
 * @return Returns 0 if success, -1 otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int start_close_timeout(int idx)
{
	if (idx < 0 || idx >= MAX_CLIENTS)
		return (-1);

	pthread_mutex_lock(&client_socks[idx].mtx_state);

	if (client_socks[idx].state != WS_STATE_OPEN)
		goto out;

	client_socks[idx].state = WS_STATE_CLOSING;

	if (pthread_create(
			&client_socks[idx].thrd_tout, NULL, close_timeout, &client_socks[idx]))
	{
		pthread_mutex_unlock(&client_socks[idx].mtx_state);
		panic("Unable to create timeout thread\n");
	}
	client_socks[idx].close_thrd = true;
out:
	pthread_mutex_unlock(&client_socks[idx].mtx_state);
	return (0);
}
#endif

/**
 * @brief Gets the IP address relative to a file descriptor opened
 * by the server.
 *
 * @param fd File descriptor target.
 *
 * @return Pointer the ip address, or NULL if fails.
 *
 * @note It is up the caller to free the returned string.
 */
char *ws_getaddress(int fd)
{
	struct sockaddr_in addr;
	socklen_t addr_size;
	char *client;

	addr_size = sizeof(struct sockaddr_in);
	if (getpeername(fd, (struct sockaddr *)&addr, &addr_size) < 0)
		return (NULL);

	client = malloc(sizeof(char) * INET_ADDRSTRLEN);
	if (!client)
		return (NULL);

	if (!inet_ntop(AF_INET, &addr.sin_addr, client, INET_ADDRSTRLEN))
	{
		free(client);
		return (NULL);
	}
	return (client);
}

#ifndef MSG_NOSIGNAL
/* dirty fix for MacOSX, will not block sigpipe, this should be fixed
   properly, maybe using something like this:
   https://gitlab.freedesktop.org/spice/spice/-/issues/39, Mac users
   are welcome to submit a patch... */
# define MSG_NOSIGNAL 0
#endif

/* Block SIGPIPE as explained here:
   https://stackoverflow.com/questions/108183/how-to-prevent-sigpipes-or-handle-them-properly */

ssize_t ws_send(int sockfd, const void *buf, size_t len)
{
  return send (sockfd, buf, len, MSG_NOSIGNAL);
}

/**
 * @brief Creates and send a WebSocket frame.
 *
 * @param fd        Target to be sent.
 * @param size      Binary message size
 * @param broadcast Enable/disable broadcast.
 * @param type      Frame type.
 *
 * @return Returns the number of bytes written, -1 if error.
 */
int ws_sendframe_init(int fd, ssize_t size, bool broadcast, int type)
{
	unsigned char frame[10]; /* Frame.             */
	uint8_t idx_first_rData; /* Index data.        */
	uint64_t length;         /* Message length.    */
	ssize_t output;          /* Bytes sent.        */
	int sock;                /* File Descript.     */

	frame[0] = (WS_FIN | type);
	length   = size;

	/* Split the size between octets. */
	if (length <= 125)
	{
		frame[1]        = length & 0x7F;
		idx_first_rData = 2;
	}

	/* Size between 126 and 65535 bytes. */
	else if (length >= 126 && length <= 65535)
	{
		frame[1]        = 126;
		frame[2]        = (length >> 8) & 255;
		frame[3]        = length & 255;
		idx_first_rData = 4;
	}

	/* More than 65535 bytes. */
	else
	{
		frame[1]        = 127;
		frame[2]        = (unsigned char)((length >> 56) & 255);
		frame[3]        = (unsigned char)((length >> 48) & 255);
		frame[4]        = (unsigned char)((length >> 40) & 255);
		frame[5]        = (unsigned char)((length >> 32) & 255);
		frame[6]        = (unsigned char)((length >> 24) & 255);
		frame[7]        = (unsigned char)((length >> 16) & 255);
		frame[8]        = (unsigned char)((length >> 8) & 255);
		frame[9]        = (unsigned char)(length & 255);
		idx_first_rData = 10;
	}

	output                 = fd >= 0 ? ws_send(CLI_SOCK(fd), frame, idx_first_rData) : 0;
	if (broadcast)
	{

		for (int i = 0; i < MAX_CLIENTS; i++)
		{
			sock = client_socks[i].client_sock;
			if ((sock > -1) && (sock != fd && sock != -fd))
				output += ws_send(CLI_SOCK(sock), frame, idx_first_rData);
		}
	}

	return ((int)output);
}

/**
 * @brief Creates and send an WebSocket frame with some payload data.
 *
 * This routine is intended to be used to create a websocket frame for
 * a given type e sending to the client. For higher level routines,
 * please check @ref ws_sendframe_txt and @ref ws_sendframe_bin.
 *
 * @param fd        Target to be sent.
 * @param msg       Message to be sent.
 * @param size      Binary message size (set it <0 for text message)
 * @param broadcast Enable/disable broadcast.
 * @param type      Frame type.
 *
 * @return Returns the number of bytes written, -1 if error.
 *
 * @note If @p size is -1, it is assumed that a text frame is being sent,
 * otherwise, a binary frame. In the later case, the @p size is used.
 */
int ws_sendframe(int fd, const char *msg, ssize_t size, bool broadcast, int type)
{

	/* Guess the size if not informed, perhaps a TXT frame. */
	if (size < 0)
		size = strlen((const char *)msg);

	if (ws_sendframe_init(fd, size, broadcast, type) < 0)
	  return -1;
	  
	ssize_t output = 0;
	
	if (fd >= 0)
	  output += ws_send(CLI_SOCK(fd), msg, size);
	
	if (broadcast)
	{

		for (int i = 0; i < MAX_CLIENTS; i++)
		{
			int sock = client_socks[i].client_sock;
			if ((sock > -1) && (sock != fd && sock != - fd))
				output += ws_send(CLI_SOCK(sock), msg, size);
		}
	}

	return ((int)output);
}

/**
 * @brief Sends a WebSocket text frame.
 *
 * @param fd         Target to be send.
 * @param msg        Text message to be send.
 * @param broadcast  Enable/disable broadcast (0-disable/anything-enable).
 *
 * @return Returns the number of bytes written, -1 if error.
 */
int ws_sendframe_txt(int fd, const char *msg, bool broadcast)
{
	return ws_sendframe(fd, msg, -1, broadcast, WS_FR_OP_TXT);
}

/**
 * @brief Sends a WebSocket binary frame.
 *
 * @param fd         Target to be send.
 * @param msg        Binary message to be send.
 * @param size       Message size (in bytes).
 * @param broadcast  Enable/disable broadcast (0-disable/anything-enable).
 *
 * @return Returns the number of bytes written, -1 if error.
 */
int ws_sendframe_bin(int fd, const char *msg, size_t size, bool broadcast)
{
	return ws_sendframe(fd, msg, size, broadcast, WS_FR_OP_BIN);
}

/**
 * @brief For a given @p fd, gets the current state for
 * the connection, or -1 if invalid.
 *
 * @param fd Client fd.
 *
 * @return Returns the connection state or -1 if
 * invalid @p fd.
 *
 * @see WS_STATE_CONNECTING
 * @see WS_STATE_OPEN
 * @see WS_STATE_CLOSING
 * @see WS_STATE_CLOSED
 */
int ws_get_state(int fd)
{
	int idx;

	if ((idx = get_client_index(fd)) == -1)
		return (-1);

	return (get_client_state(idx));
}

/**
 * @brief Close the client connection for the given @p fd
 * with normal close code (1000) and no reason string.
 *
 * @param fd Client fd.
 *
 * @return Returns 0 on success, -1 otherwise.
 *
 * @note If the client did not send a close frame in
 * TIMEOUT_MS milliseconds, the server will close the
 * connection with error code (1002).
 */
int ws_close_client(int fd)
{
	unsigned char clse_code[2];
	int cc;
	int i;

	/* Check if fd belongs to a connected client. */
	if ((i = get_client_index(fd)) == -1)
		return (-1);

	/*
	 * Instead of using do_close(), we use this to avoid using
	 * msg_ctrl buffer from wfd and avoid a race condition
	 * if this is invoked asynchronously.
	 */
	cc           = WS_CLSE_NORMAL;
	clse_code[0] = (cc >> 8);
	clse_code[1] = (cc & 0xFF);
	if (ws_sendframe(CLI_SOCK(fd), (const char *)clse_code, sizeof(char) * 2, false,
			WS_FR_OP_CLSE) < 0)
	{
		DEBUG("An error has occurred while sending closing frame!\n");
		return (-1);
	}

	/*
	 * Starts the timeout thread: if the client did not send
	 * a close frame in TIMEOUT_MS milliseconds, the server
	 * will close the connection with error code (1002).
	 */
	//	start_close_timeout(i); // fixme
	return (0);
}

/**
 * @brief Checks is a given opcode @p frame
 * belongs to a control frame or not.
 *
 * @param frame Frame opcode to be checked.
 *
 * @return Returns 1 if is a control frame, 0 otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static inline int is_control_frame(int frame)
{
	return (
		frame == WS_FR_OP_CLSE || frame == WS_FR_OP_PING || frame == WS_FR_OP_PONG);
}

/**
 * @brief Do the handshake process.
 *
 * @param wfd Websocket Frame Data.
 * @param p_index Client port index.
 *
 * @return Returns 0 if success, a negative number otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int do_handshake(struct ws_frame_data *wfd)
{
	char *response; /* Handshake response message. */
	char *p;        /* Last request line pointer.  */
	ssize_t n;      /* Read/Write bytes.           */

	/* Read the very first client message. */
	if ((n = read(wfd->sock, wfd->frm, sizeof(wfd->frm) - 1)) < 0)
		return (-1);

	/* Advance our pointers before the first next_byte(). */
	p = strstr((const char *)wfd->frm, "\r\n\r\n");
	if (p == NULL)
	{
		DEBUG("An empty line with \\r\\n was expected!\n");
		return (-1);
	}
	wfd->amt_read = n;
	wfd->cur_pos  = (size_t)((ptrdiff_t)(p - (char *)wfd->frm)) + 4;

	/* Get response. */
	if (get_handshake_response((char *)wfd->frm, &response) < 0)
	{
		DEBUG("Cannot get handshake response, request was: %s\n", wfd->frm);
		return (-1);
	}

	/* Valid request. */
	DEBUG("Handshaked, response: \n"
		  "------------------------------------\n"
		  "%s"
		  "------------------------------------\n",
		response);

	/* Send handshake. */
	if (ws_send(CLI_SOCK(wfd->sock), response, strlen(response)) < 0)
	{
		free(response);
		DEBUG("As error has occurred while handshaking!\n");
		return (-1);
	}

	/* clean up buffers. */
	free(response);
	return (0);
}

/**
 * @brief Sends a close frame, accordingly with the @p close_code
 * or the message inside @p wfd.
 *
 * @param wfd Websocket Frame Data.
 * @param close_code Websocket close code.
 *
 * @return Returns 0 if success, a negative number otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int do_close(struct ws_frame_data *wfd, int close_code)
{
	int cc; /* Close code.           */

	/* If custom close-code. */
	if (close_code != -1)
	{
		cc = close_code;
		goto custom_close;
	}

	/* If empty or have a close reason, just re-send. */
	if (wfd->frame_size == 0 || wfd->frame_size > 2)
		goto send;

	/* Parse close code and check if valid, if not, we issue an protocol error. */
	if (wfd->frame_size == 1)
		cc = wfd->msg_ctrl[0];
	else
		cc = ((int)wfd->msg_ctrl[0]) << 8 | wfd->msg_ctrl[1];

	/* Check if it's not valid, if so, we send a protocol error (1002). */
	if ((cc < 1000 || cc > 1003) && (cc < 1007 || cc > 1011) &&
		(cc < 3000 || cc > 4999))
	{
		cc = WS_CLSE_PROTERR;

	custom_close:
		wfd->msg_ctrl[0] = (cc >> 8);
		wfd->msg_ctrl[1] = (cc & 0xFF);

		if (ws_sendframe(CLI_SOCK(wfd->sock), (const char *)wfd->msg_ctrl,
				sizeof(char) * 2, false, WS_FR_OP_CLSE) < 0)
		{
			DEBUG("An error has occurred while sending closing frame!\n");
			return (-1);
		}
		return (0);
	}

	/* Send the data inside wfd->msg_ctrl. */
send:
	if (ws_sendframe(CLI_SOCK(wfd->sock), (const char *)wfd->msg_ctrl,
			wfd->frame_size, false, WS_FR_OP_CLSE) < 0)
	{
		DEBUG("An error has occurred while sending closing frame!\n");
		return (-1);
	}
	return (0);
}

/**
 * @brief Send a pong frame in response to a ping frame.
 *
 * Accordingly to the RFC, a pong frame must have the same
 * data payload as the ping frame, so we just send a
 * ordinary frame with PONG opcode.
 *
 * @param wfd Websocket frame data.
 *
 * @return Returns 0 if success and a negative number
 * otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int do_pong(struct ws_frame_data *wfd, size_t frame_size)
{
	if (ws_sendframe(CLI_SOCK(wfd->sock), (const char *)wfd->msg_ctrl, frame_size,
			false, WS_FR_OP_PONG) < 0)
	{
		wfd->error = 1;
		DEBUG("An error has occurred while ponging!\n");
		return (-1);
	}
	return (0);
}

/**
 * @brief Read a chunk of bytes and return the next byte
 * belonging to the frame.
 *
 * @param wfd Websocket Frame Data.
 *
 * @return Returns the byte read, or -1 if error.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static inline int next_byte(struct ws_frame_data *wfd)
{
	ssize_t n;

	/* If empty or full. */
	if (wfd->cur_pos == 0 || wfd->cur_pos == wfd->amt_read)
	{
		// SP: we read only one byte at a time to simplify
		// polling, this is not efficient for large messages received from the client
		if ((n = read(wfd->sock, wfd->frm, 1/*sizeof(wfd->frm)*/)) <= 0)
		{
			wfd->error = 1;
			DEBUG("An error has occurred while trying to read next byte\n");
			return (-1);
		}
		wfd->amt_read = (size_t)n;
		wfd->cur_pos  = 0;
	}
	return (wfd->frm[wfd->cur_pos++]);
}

/**
 * @brief Skips @p frame_size bytes of the current frame.
 *
 * @param wfd Websocket Frame Data.
 * @param frame_size Amount of bytes to be skipped.
 *
 * @return Returns 0 if success, a negative number
 * otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int skip_frame(struct ws_frame_data *wfd, size_t frame_size)
{
	size_t i;
	for (i = 0; i < frame_size; i++)
	{
		if (next_byte(wfd) == -1)
		{
			wfd->error = 1;
			return (-1);
		}
	}
	return (0);
}

/**
 * @brief Reads the current frame isolating data from control frames.
 * The parameters are changed in order to reflect the current state.
 *
 * @param wfd Websocket Frame Data.
 * @param opcode Frame opcode.
 * @param buf Buffer to be written.
 * @param frame_length Length of the current frame.
 * @param frame_size Total size of the frame (considering CONT frames)
 *                   read until the moment.
 * @param msg_idx Message index, reflects the current buffer pointer state.
 * @param masks Masks vector.
 * @param is_fin Is FIN frame indicator.
 *
 * @return Returns 0 if success, a negative number otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int read_frame(struct ws_frame_data *wfd,
	int opcode,
	unsigned char **buf,
	size_t *frame_length,
	size_t *frame_size,
	size_t *msg_idx,
	uint8_t *masks,
	int is_fin)
{
	unsigned char *tmp; /* Tmp message.     */
	unsigned char *msg; /* Current message. */
	int cur_byte;       /* Curr byte read.  */
	size_t i;           /* Loop index.      */

	msg = *buf;

	/* Decode masks and length for 16-bit messages. */
	if (*frame_length == 126)
		*frame_length = (((size_t)next_byte(wfd)) << 8) | next_byte(wfd);

	/* 64-bit messages. */
	else if (*frame_length == 127)
	{
		*frame_length =
			(((size_t)next_byte(wfd)) << 56) | /* frame[2]. */
			(((size_t)next_byte(wfd)) << 48) | /* frame[3]. */
			(((size_t)next_byte(wfd)) << 40) | (((size_t)next_byte(wfd)) << 32) |
			(((size_t)next_byte(wfd)) << 24) | (((size_t)next_byte(wfd)) << 16) |
			(((size_t)next_byte(wfd)) << 8) |
			(((size_t)next_byte(wfd))); /* frame[9]. */
	}

	*frame_size += *frame_length;

	/*
	 * Check frame size
	 *
	 * We need to limit the amount supported here, since if
	 * we follow strictly to the RFC, we have to allow 2^64
	 * bytes. Also keep in mind that this is still true
	 * for continuation frames.
	 */
	if (*frame_size > MAX_FRAME_LENGTH)
	{
		DEBUG("Current frame from client %d, exceeds the maximum\n"
			  "amount of bytes allowed (%zu/%d)!",
			wfd->sock, *frame_size + *frame_length, MAX_FRAME_LENGTH);

		wfd->error = 1;
		return (-1);
	}

	/* Read masks. */
	masks[0] = next_byte(wfd);
	masks[1] = next_byte(wfd);
	masks[2] = next_byte(wfd);
	masks[3] = next_byte(wfd);

	/*
	 * Abort if error.
	 *
	 * This is tricky: we may have multiples error codes from the
	 * previous next_bytes() calls, but, since we're only setting
	 * variables and flags, there is no major issue in setting
	 * them wrong _if_ we do not use their values, thing that
	 * we do here.
	 */
	if (wfd->error)
		return (-1);

	/*
	 * Allocate memory.
	 *
	 * The statement below will allocate a new chunk of memory
	 * if msg is NULL with size total_length. Otherwise, it will
	 * resize the total memory accordingly with the message index
	 * and if the current frame is a FIN frame or not, if so,
	 * increment the size by 1 to accommodate the line ending \0.
	 */
	if (*frame_length > 0)
	{
		if (!is_control_frame(opcode))
		{
			tmp = realloc(
				msg, sizeof(unsigned char) * (*msg_idx + *frame_length + is_fin));
			if (!tmp)
			{
				DEBUG("Cannot allocate memory, requested: %zu\n",
					(*msg_idx + *frame_length + is_fin));

				wfd->error = 1;
				return (-1);
			}
			msg  = tmp;
			*buf = msg;
		}

		/* Copy to the proper location. */
		for (i = 0; i < *frame_length; i++, (*msg_idx)++)
		{
			/* We were able to read? .*/
			cur_byte = next_byte(wfd);
			if (cur_byte == -1)
				return (-1);

			msg[*msg_idx] = cur_byte ^ masks[i % 4];
		}
	}

	/* If we're inside a FIN frame, lets... */
	if (is_fin && *frame_size > 0)
	{
		/* Increase memory if our FIN frame is of length 0. */
		if (!*frame_length && !is_control_frame(opcode))
		{
			tmp = realloc(msg, sizeof(unsigned char) * (*msg_idx + 1));
			if (!tmp)
			{
				DEBUG("Cannot allocate memory, requested: %zu\n", (*msg_idx + 1));

				wfd->error = 1;
				return (-1);
			}
			msg  = tmp;
			*buf = msg;
		}
		msg[*msg_idx] = '\0';
	}

	return (0);
}

/**
 * @brief Reads the next frame, whether if a TXT/BIN/CLOSE
 * of arbitrary size.
 *
 * @param wfd Websocket Frame Data.
 * @param idx Websocket connection index.
 *
 * @return Returns 0 if success, a negative number otherwise.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static int next_frame(struct ws_frame_data *wfd, int idx)
{
	unsigned char *msg_data; /* Data frame.                */
	unsigned char *msg_ctrl; /* Control frame.             */
	uint8_t masks_data[4];   /* Masks data frame array.    */
	uint8_t masks_ctrl[4];   /* Masks control frame array. */
	size_t msg_idx_data;     /* Current msg index.         */
	size_t msg_idx_ctrl;     /* Current msg index.         */
	size_t frame_length;     /* Frame length.              */
	size_t frame_size;       /* Current frame size.        */
	uint8_t opcode;          /* Frame opcode.              */
	uint8_t is_fin;          /* Is FIN frame flag.         */
	uint8_t mask;            /* Mask.                      */
	int cur_byte;            /* Current frame byte.        */

	msg_data        = NULL;
	msg_ctrl        = wfd->msg_ctrl;
	is_fin          = 0;
	frame_length    = 0;
	frame_size      = 0;
	msg_idx_data    = 0;
	msg_idx_ctrl    = 0;
	wfd->frame_size = 0;
	wfd->frame_type = -1;
	wfd->msg        = NULL;
		
	/* Read until find a FIN or a unsupported frame. */
	do
	{
		/*
		 * Obs: next_byte() can return error if not possible to read the
		 * next frame byte, in this case, we return an error.
		 *
		 * However, please note that this check is only made here and in
		 * the subsequent next_bytes() calls this also may occur too.
		 * wsServer is assuming that the client only create right
		 * frames and we will do not have disconnections while reading
		 * the frame but just when waiting for a frame.
		 */

		cur_byte = next_byte(wfd);
		if (cur_byte == -1)
			return (-1);

		is_fin = (cur_byte & 0xFF) >> WS_FIN_SHIFT;
		opcode = (cur_byte & 0xF);

		/*
		 * Check for RSV field.
		 *
		 * Since wsServer do not negotiate extensions if we receive
		 * a RSV field, we must drop the connection.
		 */
		if (cur_byte & 0x70)
		{
			DEBUG("RSV is set while wsServer do not negotiate extensions!\n");
			wfd->error = 1;
			break;
		}

		/*
		 * Check if the current opcode makes sense:
		 * a) If we're inside a cont frame but no previous data frame
		 *
		 * b) If we're handling a data-frame and receive another data
		 *    frame. (it's expected to receive only CONT or control
		 *    frames).
		 *
		 * It is worth to note that in a), we do not need to check
		 * if the previous frame was FIN or not: if was FIN, an
		 * on_message event was triggered and this function returned;
		 * so the only possibility here is a previous non-FIN data
		 * frame, ;-).
		 */
		if ((wfd->frame_type == -1 && opcode == WS_FR_OP_CONT) ||
			(wfd->frame_type != -1 && !is_control_frame(opcode) &&
				opcode != WS_FR_OP_CONT))
		{
			DEBUG("Unexpected frame was received!, opcode: %d, previous: %d\n",
				opcode, wfd->frame_type);
			wfd->error = 1;
			break;
		}

		/* Check if one of the valid opcodes. */
		if (opcode == WS_FR_OP_TXT || opcode == WS_FR_OP_BIN ||
			opcode == WS_FR_OP_CONT || opcode == WS_FR_OP_PING ||
			opcode == WS_FR_OP_PONG || opcode == WS_FR_OP_CLSE)
		{
			/*
			 * Check our current state: if CLOSING, we only accept close
			 * frames.
			 *
			 * Since the server may, at any time, asynchronously, asks
			 * to close the client connection, we should terminate
			 * immediately.
			 */
			if (get_client_state(idx) == WS_STATE_CLOSING && opcode != WS_FR_OP_CLSE)
			{
				DEBUG(
					"Unexpected frame received, expected CLOSE (%d), received: (%d)",
					WS_FR_OP_CLSE, opcode);
				wfd->error = 1;
				break;
			}

			/* Only change frame type if not a CONT frame. */
			if (opcode != WS_FR_OP_CONT && !is_control_frame(opcode))
				wfd->frame_type = opcode;

			mask         = next_byte(wfd);
			frame_length = mask & 0x7F;
			frame_size   = 0;
			msg_idx_ctrl = 0;

			/*
			 * We should deny non-FIN control frames or that have
			 * more than 125 octets.
			 */
			if (is_control_frame(opcode) && (!is_fin || frame_length > 125))
			{
				DEBUG("Control frame bigger than 125 octets or not a FIN frame!\n");
				wfd->error = 1;
				break;
			}

			/* Normal data frames. */
			if (opcode == WS_FR_OP_TXT || opcode == WS_FR_OP_BIN ||
				opcode == WS_FR_OP_CONT)
			{
				read_frame(wfd, opcode, &msg_data, &frame_length, &wfd->frame_size,
					&msg_idx_data, masks_data, is_fin);
			}

			/*
			 * We _never_ send a PING frame, so it's not expected to receive a PONG
			 * frame. However, the specs states that a client could send an
			 * unsolicited PONG frame. The server just have to ignore the
			 * frame.
			 *
			 * The skip amount will always be 4 (masks vector size) + frame size
			 */
			else if (opcode == WS_FR_OP_PONG)
			{
				skip_frame(wfd, 4 + frame_length);
				is_fin = 0;
				continue;
			}

			/* We should answer to a PING frame as soon as possible. */
			else if (opcode == WS_FR_OP_PING)
			{
				if (read_frame(wfd, opcode, &msg_ctrl, &frame_length, &frame_size,
						&msg_idx_ctrl, masks_ctrl, is_fin) < 0)
					break;

				if (do_pong(wfd, frame_size) < 0)
					break;

				/* Quick hack to keep our loop. */
				is_fin = 0;
			}

			/* We interrupt the loop as soon as we find a CLOSE frame. */
			else
			{
				if (read_frame(wfd, opcode, &msg_ctrl, &frame_length, &frame_size,
						&msg_idx_ctrl, masks_ctrl, is_fin) < 0)
					break;

				/* Since we're aborting, we can scratch the 'data'-related
				 * vars here. */
				wfd->frame_size = frame_size;
				wfd->frame_type = WS_FR_OP_CLSE;
				free(msg_data);
				return (0);
			}
		}

		/* Anything else (unsupported frames). */
		else
		{
			DEBUG("Unsupported frame opcode: %d\n", opcode);
			/* We should consider as error receive an unknown frame. */
			wfd->frame_type = opcode;
			wfd->error      = 1;
		}

	} while (!is_fin && !wfd->error);

	/* Check for error. */
	if (wfd->error)
	{
		free(msg_data);
		wfd->msg = NULL;
		return (-1);
	}

	wfd->msg = msg_data;
	return (0);
}

static void close_connection(int connection_index)
{
	if (client_socks[connection_index].state != WS_STATE_CLOSED)
	{
		/* Removes client socket from socks list. */
		close(client_socks[connection_index].client_sock);
		client_socks[connection_index].client_sock = -1;
		client_socks[connection_index].state       = WS_STATE_CLOSED;
	}
}

static struct ws_message * append_message (struct ws_message * messages,
					   struct ws_message msg)
{
	if (messages == NULL) {
		messages = malloc (sizeof(struct ws_message));
		messages[0].fd = -1;
	}
	int n = 0;
	while (messages[n].fd >= 0) n++;
	messages = realloc (messages, (n + 2)*sizeof(struct ws_message));
	messages[n] = msg;
	messages[n + 1].fd = -1;
	return messages;
}

static struct ws_message * handle_message (int connection_index, struct ws_message * messages)
{
	int sock = client_socks[connection_index].client_sock;
	struct ws_frame_data * wfd = &client_socks[connection_index].frame_data;
	
	/* Read next frame until client disconnects or an error occur. */
	int status = next_frame(wfd, connection_index);
	if (status >= 0)
	{
		/* Text/binary event. */
		if ((wfd->frame_type == WS_FR_OP_TXT || wfd->frame_type == WS_FR_OP_BIN) &&
		    !wfd->error)
		{
			messages = append_message (messages, (struct ws_message) {
					sock, (char *) wfd->msg, wfd->frame_size, wfd->frame_type });
		}
				
		/* Close event. */
		else if (wfd->frame_type == WS_FR_OP_CLSE && !wfd->error)
		{
			status = -1;
					
			/*
			 * We only send a CLOSE frame once, if we're already
			 * in CLOSING state, there is no need to send.
			 */
			if (get_client_state(connection_index) != WS_STATE_CLOSING)
			{
				set_client_state(connection_index, WS_STATE_CLOSING);
						
				/* We only send a close frameSend close frame */
				do_close(wfd, -1);
			}

			free(wfd->msg);
		}

	}
		
	if (status < 0) {

		messages = append_message (messages, (struct ws_message) {sock, NULL, 0, WS_FR_OP_CLSE });
	
		close_connection(connection_index);
	}

	return messages;
}

/**
 * @brief Establishes to connection with the client and trigger
 * events when occurs one.
 *
 * @param vsock Client connection index.
 * @param messages The array of messages.
 *
 * @return Returns the array of connection messages.
 *
 * @attention This is part of the internal API and is documented just
 * for completeness.
 */
static struct ws_message * ws_establishconnection(int connection_index, struct ws_message * messages)
{
	struct ws_frame_data * wfd; /* WebSocket frame data.   */
	int sock;                 /* File descriptor.        */

	sock             = client_socks[connection_index].client_sock;
	wfd              = &client_socks[connection_index].frame_data;
	
	/* Prepare frame data. */
	memset(wfd, 0, sizeof(struct ws_frame_data));
	wfd->sock = sock;

	/* Do handshake. */
	if (do_handshake(wfd) < 0) {
		close_connection(connection_index);
		return messages;
	}

	/* Change state. */
	set_client_state(connection_index, WS_STATE_OPEN);

	return append_message (messages, (struct ws_message) {sock, NULL, 0, WS_FR_OP_OPEN});
}

/**
 * @brief Opens a WebSocket.
 *
 * @param port Server port.
 *
 * @return the socket.
 */
int ws_socket_open (uint16_t port)
{
	int sock;                  /* Current socket.        */
	struct sockaddr_in server; /* Server.                */

	/* Create socket. */
	sock = socket(AF_INET, SOCK_STREAM, 0);
	if (sock < 0)
		panic("Could not create socket");

	/* Reuse previous address. */
	if (setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &(int){1}, sizeof(int)) < 0)
		panic("setsockopt(SO_REUSEADDR) failed");

	/* Prepare the sockaddr_in structure. */
	server.sin_family      = AF_INET;
	server.sin_addr.s_addr = INADDR_ANY;
	server.sin_port        = htons(port);

	/* Bind. */
	if (bind(sock, (struct sockaddr *)&server, sizeof(server)) < 0) {
		close(sock);
		return -1;
	}

	/* Listen. */
	listen(sock, MAX_CLIENTS);

	memset(client_socks, -1, sizeof(client_socks));

	return sock;
}

/**
 * @brief Polls a WebSocket.
 *
 * @param sock The WebSocket (as returned by ws_socket_open()).
 * @param timeout The poll timeout (in milliseconds).
 *
 * @return the array of messages or NULL.
 */
struct ws_message * ws_socket_poll (int sock, int timeout)
{
	int nr = 0;
	struct ws_message * messages = NULL;

	do {

		struct pollfd fds[MAX_CLIENTS + 1] = {{sock, POLLIN, 0}};
		for (int i = 0; i < MAX_CLIENTS; i++)
		{
			fds[i + 1].fd = client_socks[i].client_sock;
			fds[i + 1].events = POLLIN;
		}

		if (poll (fds, MAX_CLIENTS + 1, timeout) < 0)
			panic("Error on polling clients");
		timeout = 0;

		nr = 0;
		for (int j = 0; j < MAX_CLIENTS; j++)
		{

			if (fds[j].revents == POLLIN)
			{

				nr++;
			
				if (j == 0) /* Accept */
				{

					/* Adds client socket to socks list. */
					int connection_index = -1;
					for (int i = 0; i < MAX_CLIENTS; i++)
					{
						if (client_socks[i].client_sock == -1)
						{
							connection_index            = i;
							break;
						}
					}

					if (connection_index >= 0) {
					
						struct sockaddr_in client; /* Client.                */
						int len = sizeof(struct sockaddr_in);
						int new_sock = accept(sock, (struct sockaddr *)&client, (socklen_t *)&len);
						if (new_sock < 0)
							panic("Error on accepting connections..");
					
						client_socks[connection_index].client_sock = new_sock;
						client_socks[connection_index].state       = WS_STATE_CONNECTING;

						messages = ws_establishconnection (connection_index, messages);
					}
				}
				else /* Message */
				{
					
					messages = handle_message (j - 1, messages);
					
				}
				
			}
		}
		
	} while (nr);
	
	return messages;
}


#ifdef AFL_FUZZ
/**
 * @brief WebSocket fuzzy test routine
 *
 * @param evs  Events structure.
 *
 * @param file File to be read.
 *
 * @return Returns 0, or crash.
 *
 * @note This is a 'fuzzy version' of the function @ref ws_socket.
 * This routine do not listen to any port nor accept multiples
 * connections. It is intended to read a stream of frames through a
 * file and process it as if they are coming from a socket.
 *
 * This behavior enables us to test wsServer against fuzzers, like
 * AFL, and see if it crashes, hangs or behaves normally, even under
 * weird conditions.
 */
int ws_file(struct ws_events *evs, const char *file)
{
	int sock;
	sock = open(file, O_RDONLY);
	if (sock < 0)
		panic("Invalid file\n");

	/* Copy events. */
	memcpy(&ports[0].events, evs, sizeof(struct ws_events));
	ports[0].port_number = 0;

	/* Clear client socks list. */
	memset(client_socks, -1, sizeof(client_socks));

	/* Set client settings. */
	client_socks[0].client_sock = sock;
	client_socks[0].port_index  = 0;

	ws_establishconnection((void *)(intptr_t)0);
	return (0);
}
#endif
