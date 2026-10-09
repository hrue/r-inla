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
		switch (color) {
		case COLOR_STRONG_RED:
			fprintf(stderr, "\x1b[1;31m");
			break;
		case COLOR_STRONG_BLUE:
			fprintf(stderr, "\x1b[94m");
			break;
		case COLOR_RESET:
		default:
			fprintf(stderr, "\x1b[0m");
			break;
		}
	}
}

// 3. The Core Quality Logging Macros
// Strong Red Logging Wrapper
#define LOG_ERROR(fmt, ...) do {					\
		set_terminal_color(COLOR_STRONG_RED);			\
		fprintf(stderr, "***ERROR*** " fmt "\n", ##__VA_ARGS__); \
		set_terminal_color(COLOR_RESET);			\
	} while(0)

// Strong Blue Logging Wrapper
#define LOG_WARNING(fmt, ...) do {						\
		set_terminal_color(COLOR_STRONG_BLUE);			\
		fprintf(stderr, "***INFO***  " fmt "\n", ##__VA_ARGS__); \
		set_terminal_color(COLOR_RESET);			\
	} while(0)

int main()
{
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
