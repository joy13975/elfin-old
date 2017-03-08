To compile:
	make

To run:
	./elfin <arguments to override settings.json>

To edit settings:
	<your editor> settings.json

Notes:
	1. Makefile adds jutil/src/ to include path, that's why util.h is included as
	if it was in the same directory everywhere