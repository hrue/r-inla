#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <wchar.h>

#define SIZE 10

int main(void)
{
	setlocale(LC_ALL, "");
	wchar_t buf[SIZE+1];
	wchar_t *pat = L"привет мир";
	wchar_t str[SIZE+2];

	FILE *f1;
	FILE *f2;

	f1 = fopen("./вход","r");
	f2 = fopen("./выход","w");

	fgetws(buf, SIZE+1, f1);

	if (wcsncmp(buf, pat, SIZE) == 0) {
		swprintf(str, SIZE+2, L"% 11ls", buf);
		fputws(str, f2);
	}

	fclose(f1);
	fclose(f2);

	exit(0);
}
