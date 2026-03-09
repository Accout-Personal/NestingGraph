#pragma once

// Install handlers so Ctrl+C / terminal close exits with `exit_code` (default 130).
void install_cancel_exit_handlers(int exit_code = 130);