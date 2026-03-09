#include "cancel_exit.hpp"

#ifdef _WIN32
  #include <windows.h>

  static int g_exit_code = 130;

  static BOOL WINAPI ConsoleHandler(DWORD type) {
    switch (type) {
      case CTRL_C_EVENT:
      case CTRL_BREAK_EVENT:
      case CTRL_CLOSE_EVENT:   // user clicks X on console
      case CTRL_LOGOFF_EVENT:
      case CTRL_SHUTDOWN_EVENT:
        ExitProcess((UINT)g_exit_code); // immediate, no cleanup
        return TRUE;
      default:
        return FALSE;
    }
  }

  void install_cancel_exit_handlers(int exit_code) {
    g_exit_code = exit_code;
    SetConsoleCtrlHandler(ConsoleHandler, TRUE);
  }

#else
  #include <signal.h>
  #include <unistd.h>

  static int g_exit_code = 130;

  static void on_signal(int) {
    _exit(g_exit_code); // immediate, async-signal-safe
  }

  void install_cancel_exit_handlers(int exit_code) {
    g_exit_code = exit_code;

    struct sigaction sa{};
    sa.sa_handler = on_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;

    sigaction(SIGINT,  &sa, nullptr); // Ctrl+C
    sigaction(SIGHUP,  &sa, nullptr); // terminal closed
    sigaction(SIGTERM, &sa, nullptr); // kill/terminate
  }
#endif