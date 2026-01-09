#ifndef _TEDSBELLHOP_EXPORT_H_
#define _TEDSBELLHOP_EXPORT_H_



#if defined(_WIN32)
#if defined(TEDSBELLHOP_EXPORTS)
#define TEDSBELLHOP_API __declspec(dllexport)
#else
#define TEDSBELLHOP_API __declspec(dllimport)
#endif
#else
#define TEDSBELLHOP_API __attribute__((visibility("default")))
#endif



#endif // !_TEDSBELLHOP_EXPORT_H_
