#include "GMRFLib/GMRFLib.h"

#if defined(_WIN32) || defined(_WIN64)
#       include <windows.h>
#       define PLATFORM_WINDOWS 1
#else
#       define PLATFORM_WINDOWS 0
#endif

#if defined(INLA_WITH_COLOR)

void GMRFLib_color_core(GMRFLib_color_tp color)
{
	GMRFLib_color_set(color);
	fprintf(stderr, "\n\n");
	fflush(stderr);
}

void GMRFLib_color_start_error(void)
{
	GMRFLib_color_core(COLOR_BG_RED_FG_WHITE);
}

void GMRFLib_color_start_warning(void)
{
	GMRFLib_color_core(COLOR_BG_BLUE_FG_WHITE);
}

void GMRFLib_color_start_info(void)
{
	GMRFLib_color_core(COLOR_BG_GREEN_FG_WHITE);
}

void GMRFLib_color_reset(void)
{
	GMRFLib_color_core(COLOR_RESET);
}

int windows_ansi_ready = 0;
void GMRFLib_color_init(void)
{
#       if PLATFORM_WINDOWS
	int out_ok = 0;
	int err_ok = 0;

	// Enable ANSI processing for Standard Output (stdout)
	HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
	if (hOut != INVALID_HANDLE_VALUE) {
		DWORD dwMode = 0;
		if (GetConsoleMode(hOut, &dwMode)) {
			dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
			if (SetConsoleMode(hOut, dwMode))
				out_ok = 1;
		}
	}
	// Enable ANSI processing for Standard Error (stderr)
	HANDLE hErr = GetStdHandle(STD_ERROR_HANDLE);
	if (hErr != INVALID_HANDLE_VALUE) {
		DWORD dwMode = 0;
		if (GetConsoleMode(hErr, &dwMode)) {
			dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
			if (SetConsoleMode(hErr, dwMode))
				err_ok = 1;
		}
	}
	// If both handles support ANSI escape codes, we can use uniform string sequences!
	if (out_ok && err_ok) {
		windows_ansi_ready = 1;
	}
#       endif
}

void GMRFLib_color_set(GMRFLib_color_tp color)
{
	// If we are on Windows and modern ANSI setup failed, fall back to legacy console attributes
	if (PLATFORM_WINDOWS && !windows_ansi_ready) {
#       if defined(_WIN32) || defined(_WIN64)
		HANDLE hConsole = GetStdHandle(STD_ERROR_HANDLE);
		switch (color) {
		case COLOR_BG_RED_FG_WHITE:
			SetConsoleTextAttribute(hConsole,
						BACKGROUND_RED | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
			break;
		case COLOR_BG_BLUE_FG_WHITE:
			SetConsoleTextAttribute(hConsole,
						BACKGROUND_BLUE | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
			break;
		case COLOR_BG_GREEN_FG_WHITE:
			SetConsoleTextAttribute(hConsole,
						BACKGROUND_GREEN | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
			break;
		case COLOR_RESET:
		default:
			SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
			break;
		}
#       endif
	} else {
		// Use modern high-intensity ANSI strings (Works on macOS, Linux, and Modern Windows CMD)
		switch (color) {
		case COLOR_BG_RED_FG_WHITE:
			// ID 124 triggers an absolute dark crimson/brick red block
			fprintf(stderr, "\x1b[48;5;124m\x1b[1;97m");
			break;
		case COLOR_BG_GREEN_FG_WHITE:
			// ID 22 forces an absolute dark forest green block
			fprintf(stderr, "\x1b[48;5;22m\x1b[1;97m");
			break;
		case COLOR_BG_BLUE_FG_WHITE:
			// ID 19 forces a beautiful deep navy blue block
			fprintf(stderr, "\x1b[48;5;19m\x1b[1;97m");
			break;
		case COLOR_RESET:
		default:
			fprintf(stderr, "\x1b[0m");	       // Full style reset
			break;
		}
	}
}

#else

//  no colors, just define empty functions

void GMRFLib_color_start_error(void)
{
}

void GMRFLib_color_start_warning(void)
{
}

void GMRFLib_color_start_info(void)
{
}

void GMRFLib_color_reset(void)
{
}

void GMRFLib_color_init(void)
{
}

void GMRFLib_color_set(GMRFLib_color_tp UNUSED(color))
{
}

#endif
