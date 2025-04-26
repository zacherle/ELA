#ifndef VERSION_H
#define VERSION_H

#define VERS "v0.10.3"

// Version information header file
#ifdef __cplusplus
#include <iostream>
#include <string>

std::string getVersion() {
    return (VERS);
}
std::string getBuildDate() {
    return __DATE__;
}
std::string getBuildTime() {
    return __TIME__;
}
std::string getBuildInfo() {
    return "ELA " + getVersion() + " | Date: " + getBuildDate() + " | Time: " + getBuildTime();
}
/*
std::string getBuildInfoWithCommit() {
    return getBuildInfo() + " | Commit: " + std::string(__GIT_COMMIT__);
}
std::string getBuildInfoWithBranch() {
    return getBuildInfo() + " | Branch: " + std::string(__GIT_BRANCH__);
}
std::string getBuildInfoWithCommitAndBranch() {
    return getBuildInfo() + " | Commit: " + std::string(__GIT_COMMIT__) + " | Branch: " + std::string(__GIT_BRANCH__);
}
*/
#else

//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
//#include <time.h>
char* getVersion() {
    return (VERS);
}
char* getBuildDate() {
    return __DATE__;
}
char* getBuildTime() {
    return __TIME__;
}
char* getBuildInfo() {
    static char buildInfo[256];
    snprintf(buildInfo, sizeof(buildInfo), "ELA %s %s %s", getVersion(), getBuildDate(), getBuildTime());
    return buildInfo;
}

#endif

#endif // VERSION_H
