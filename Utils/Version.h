#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 ""
#define GIT_REFSPEC "refs/heads/master"
#define GIT_LOCAL_STATUS "HEAD-HASH-NOTFOUND"

#ifdef DL_OUTPUT

#endif

#endif
