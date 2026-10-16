#include <stdio.h>

#if defined(_WIN32) || defined(_WIN64)
#       include <windows.h>
#       define PLATFORM_WINDOWS 1
#else
#       define PLATFORM_WINDOWS 0
#endif

// 1. Text Styling Enumeration
typedef enum {
	COLOR_RESET,
	COLOR_STRONG_RED,
	COLOR_STRONG_BLUE
} TerminalColor;

// 2. Underlying Color Switcher Function
void set_terminal_color(TerminalColor color)
{
	if (PLATFORM_WINDOWS) {
#if defined(_WIN32) || defined(_WIN64)
		HANDLE hConsole = GetStdHandle(STD_ERROR_HANDLE);
		switch (color) {
		case COLOR_STRONG_RED:
			SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_INTENSITY);
			break;
		case COLOR_STRONG_BLUE:
			SetConsoleTextAttribute(hConsole, FOREGROUND_BLUE | FOREGROUND_INTENSITY);
			break;
		case COLOR_RESET:
		default:
			SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
			break;
		}
#endif
	} else {
		// Use modern high-intensity ANSI strings (Works on macOS, Linux, and Modern Windows CMD)
		switch (color) {
		case COLOR_BG_RED_FG_WHITE:
			fprintf(stderr, "\x1b[41m\x1b[1;97m"); // Red Background + Bold White Text
			break;
		case COLOR_BG_BLUE_FG_WHITE:
			fprintf(stderr, "\x1b[44m\x1b[1;97m"); // Blue Background + Bold White Text
			break;
		case COLOR_RESET:
		default:
			fprintf(stderr, "\x1b[0m");           // Full style reset
			break;
		}
	}
}

// 3. The Core Quality Logging Macros
// Error Macro: Full red background block with bold white text
#define LOG_ERROR(fmt, ...) do {					\
		set_terminal_color(COLOR_BG_RED_FG_WHITE);		\
		fprintf(stderr, " ERROR  " fmt " ", ##__VA_ARGS__);	\
		set_terminal_color(COLOR_RESET);			\
		fprintf(stderr, "\n");					\
	} while(0)

// Warning/Info Macro: Full blue background block with bold white text
#define LOG_WARNING(fmt, ...) do {					\
		set_terminal_color(COLOR_BG_BLUE_FG_WHITE);		\
		fprintf(stderr, "  INFO  " fmt " ", ##__VA_ARGS__);	\
		set_terminal_color(COLOR_RESET);			\
		fprintf(stderr, "\n");					\
	} while(0)

int main()
{
	// Initialize Windows Console Mode for ANSI sequences if applicable
	if (PLATFORM_WINDOWS) {
		init_terminal_processing();
	}

	// Dynamic initialization variables for simulation
	char *subsystem = "GMRFLib";
	int connection_id = 1047;
	double execution_time = 0.042;

	// --- Test Phase 1: Informational Logging ---
	LOG_WARNING("Initializing %s subsystem...", subsystem);
	LOG_WARNING("Connection established safely on channel #%d.", connection_id);

	// --- Test Phase 2: Critical Error Logging ---
	LOG_ERROR("Process aborted after %.3f seconds due to memory bounds breach.", execution_time);
	LOG_ERROR("Unable to write matrix values to target pointer memory map.");

	// --- Test Phase 3: Verification Output ---
	printf("Standard system validation routine complete. Color tracking restored.\n");

	return 0;
}
