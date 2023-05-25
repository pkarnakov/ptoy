#pragma once

#define PTOY_CAT(x, y) x##y
#define PTOY_XCAT(x, y) PTOY_CAT(x, y)
#define USEFLAG(x) PTOY_XCAT(0, USE_##x)
