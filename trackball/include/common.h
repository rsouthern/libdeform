#ifndef _common_h_CAL_
#define _common_h_CAL_

#include <iostream>

#define CLAMP(val,min,max) (((val)<(min))?(min):((val)>(max))?(max):(val))
#define disp_err(msg) cout << __FILE__ << ':' << __LINE__ << ": " << msg << endl;

#endif
