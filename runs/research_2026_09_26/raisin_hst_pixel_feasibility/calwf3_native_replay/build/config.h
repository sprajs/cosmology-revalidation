#ifndef HSTCAL_CONFIG_H
#define HSTCAL_CONFIG_H
#define HAVE_LIMITS__PATH_MAX
/* #undef HAVE_SYS_LIMITS__PATH_MAX */
/* #undef HAVE_SYS_SYSLIMITS__PATH_MAX */

#if defined(HAVE_LIMITS__PATH_MAX)
#include <limits.h>
#elif defined(HAVE_SYS_LIMITS__PATH_MAX)
#include <sys/limits.h>
#elif defined(HAVE_SYS_SYSLIMITS__PATH_MAX)
#include <sys/syslimits.h>
#endif

#define HSTCAL_VERSION "unknown"
#define HSTCAL_VERSION_BRANCH "unknown"
#define HSTCAL_VERSION_COMMIT "unknown"

#endif // HSTCAL_CONFIG_H
