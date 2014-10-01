#!/usr/bin/python

''' Class mycolor '''

class mycolor:
    ESCAPE = '\033[%sm' # escape
    BOLD = ESCAPE % '1;' # bold 
    UNDERLINE = ESCAPE % '4;' # underline
    # COLORS
    HEADER = '\033[95m' # pink
    BLUE = '\033[94m' # blue
    GREEN = '\033[92m' # green
    WARNING = '\033[93m' # yellow
    FAIL = '\033[91m' # red
    ENDC = ESCAPE % '0' # end escape sequence
