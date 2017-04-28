#include <sys/resource.h>

void setstacksizeunlimited(void)
{
   struct rlimit limit;
   getrlimit(RLIMIT_STACK, &limit);

   limit.rlim_cur=limit.rlim_max;

   setrlimit(RLIMIT_STACK, &limit);
}
