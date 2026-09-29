// The application instantiates NanoSVG in BitmapCache.cpp. Supply the same
// implementation to this headless exporter test without linking the GUI.
#define NANOSVG_IMPLEMENTATION
#include "nanosvg/nanosvg.h"
