#ifndef VERSION_H
#define VERSION_H

// Version information header file
#ifdef __cplusplus
extern "C" {
#endif
char* getVersion();
char* getBuildDate();
char* getBuildTime();
char* getBuildInfo();
#ifdef __cplusplus
}
#endif
#endif // VERSION_H
