#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
    #include <winsock2.h>
    #include <ws2tcpip.h>
    #pragma comment(lib, "ws2_32.lib")
    #define CLOSE_SOCKET(s) closesocket(s)
    #define CLEANUP_SOCKETS() WSACleanup()
#else
    #include <unistd.h>
    #include <arpa/inet.h>
    #include <sys/socket.socket.h>
    #define CLOSE_SOCKET(s) close(s)
    #define CLEANUP_SOCKETS()
#endif

#define SERVER_IP "79.177.133.117"  // Replace with your actual Public IP
#define PORT 65432

int main(void) {
#ifdef _WIN32
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
        printf("WSAStartup failed.\n");
        return 1;
    }
#endif

    int sock = 0;
    struct sockaddr_in serv_addr;
    char *message = "Hello World";

    // 1. Create socket
    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        printf("Socket creation error\n");
        CLEANUP_SOCKETS();
        return -1;
    }

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(PORT);

    // 2. Convert IP address from text to binary form
    if (inet_pton(AF_INET, SERVER_IP, &serv_addr.sin_addr) <= 0) {
        printf("Invalid address / Address not supported\n");
        CLOSE_SOCKET(sock);
        CLEANUP_SOCKETS();
        return -1;
    }

    // 3. Connect to the Python server
    printf("Connecting to %s:%d...\n", SERVER_IP, PORT);
    if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
        printf("Connection Failed. Check IP, Port Forwarding, or Firewall.\n");
        CLOSE_SOCKET(sock);
        CLEANUP_SOCKETS();
        return -1;
    }

    // 4. Send message
    send(sock, message, strlen(message), 0);
    printf("Message sent: %s\n", message);

    // 5. Clean up
    CLOSE_SOCKET(sock);
    CLEANUP_SOCKETS();

    return 0;
}