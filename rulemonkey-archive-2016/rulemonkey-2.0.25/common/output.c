/*
** Copyright 2008 Northern Arizona University and
** The Translational Genomics Research Institute
**
**  This file is part of Rulemonkey.
**
**  Rulemonkey is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Rulemonkey is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Rulemonkey.  If not, see <http://www.gnu.org/licenses/>.
**
*/


#include "stdio.h"
#include "stdarg.h"
#include "output.h"
#include "constants.h"
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

void suppress_warnings4() {SUPPRESS_WARNINGS;}

// Global variables.
int current_debug_level;
FILE *my_out;
FILE *my_bng_out;
FILE *my_err;
FILE *my_dbg;

// Local static variables.
static int my_out_open;
static int my_bng_out_open;
static int my_err_open;
static int my_dbg_open;

#define OUT_SUFFIX ".out"
#define BNG_OUT_SUFFIX ".gdat"
#define ERR_SUFFIX ".error"
#define DBG_SUFFIX ".info"

int
output_begin(const char *name)
{
  char *raw_name;
  char *pos;
  size_t raw_name_size;
  size_t bng_out_name_size;
  size_t out_name_size;
  size_t err_name_size;
  size_t dbg_name_size;
  char *bng_out_name;
  char *out_name;
  char *err_name;
  char *dbg_name;

  if (NULL == name)
  {
    return 0;
  }

  current_debug_level = DEFAULT_DEBUG_LEVEL;

  raw_name = strdup(name);
  pos = strrchr(raw_name, '.');
  if (pos != NULL)
  {
    // Make sure we will not overwrite the file with our output files.
    if ((0 == strcmp(pos, OUT_SUFFIX)) ||
        (0 == strcmp(pos, BNG_OUT_SUFFIX)) ||
        (0 == strcmp(pos, ERR_SUFFIX)) ||
        (0 == strcmp(pos, DBG_SUFFIX)))
    {
      free(raw_name);

      fprintf(stderr, "Filename suffix '%s' conflicts with output suffix\n",
              pos);

      return 0;
    }

    // Remove extension from name.
    *pos = '\0';
  } 


  raw_name_size = strlen(raw_name);

  // sizeof() already includes NUL character at end of string.
  out_name_size = raw_name_size + sizeof(OUT_SUFFIX);
  bng_out_name_size = raw_name_size + sizeof(BNG_OUT_SUFFIX);
  err_name_size = raw_name_size + sizeof(ERR_SUFFIX);
  dbg_name_size = raw_name_size + sizeof(DBG_SUFFIX);
  out_name = calloc(sizeof(char), out_name_size);
  bng_out_name = calloc(sizeof(char), bng_out_name_size);
  err_name = calloc(sizeof(char), err_name_size);
  dbg_name = calloc(sizeof(char), dbg_name_size);

  strcpy(out_name, raw_name);
  strcpy(bng_out_name, raw_name);
  strcpy(err_name, raw_name);
  strcpy(dbg_name, raw_name);

  strcat(out_name, OUT_SUFFIX);
  strcat(bng_out_name, BNG_OUT_SUFFIX);
  strcat(err_name, ERR_SUFFIX);
  strcat(dbg_name, DBG_SUFFIX);


  my_out = fopen(out_name, "w");
  my_out_open = 1;
  my_bng_out = fopen(bng_out_name, "w");
  my_bng_out_open = 1;

  /* Write errors to stderr. */
  my_err = stderr;
  my_err_open = 0;
  //my_err = fopen(err_name, "w");

  my_dbg = fopen(dbg_name, "w");
  my_dbg_open = 1;

#ifdef OUTPUT_FLUSH
  setbuf(my_out, NULL);
  setbuf(my_bng_out, NULL);
  setbuf(my_err, NULL);
  setbuf(my_dbg, NULL);
#endif

  free(raw_name);
  free(out_name);
  free(bng_out_name);
  free(err_name);
  free(dbg_name);

  return 1;
}

void output_end()
{
  if (1 == my_out_open) fclose(my_out);
  if (1 == my_bng_out_open) fclose(my_out);
  if (1 == my_err_open) fclose(my_err);
  if (1 == my_dbg_open) fclose(my_dbg);
}

void
set_debug_level(const int debug_level)
{
  current_debug_level = debug_level;
}

FILE *
custom_get_fh(const output_type ptype)
{
  switch(ptype)
  {
    case OUTPUT_STANDARD: 
      return my_out;
      break;
    case OUTPUT_BNG_STANDARD: 
      return my_bng_out;
      break;
    case OUTPUT_DEBUG: 
      return my_dbg;
      break;
    case OUTPUT_ERROR: 
      return my_err;
      break;
  }

  return stdout;
}

void
custom_printf(const output_type ptype, const int level, const char *fmt, ...)
{
  va_list ap;

  if (level > current_debug_level)
  {
    // Don't print out statements that are more verbose than desired.
    return;
  }

  va_start(ap, fmt);
  switch(ptype)
  {
    case OUTPUT_STANDARD: 
      vfprintf(my_out, fmt, ap);
      break;
    case OUTPUT_BNG_STANDARD: 
      vfprintf(my_bng_out, fmt, ap);
      break;
    case OUTPUT_DEBUG: 
      vfprintf(my_dbg, fmt, ap);
      break;
    case OUTPUT_ERROR: 
      vfprintf(my_err, fmt, ap);
      break;
  }
  va_end(ap);
}

void
bng_standard_printf(const char *fmt, ...)
{
  va_list ap;

  // Output to standard output file.
  va_start(ap, fmt);
  vfprintf(my_bng_out, fmt, ap);
  va_end(ap);
}

void
standard_printf(const char *fmt, ...)
{
  va_list ap;

  // Output to standard output file.
  va_start(ap, fmt);
  vfprintf(my_out, fmt, ap);
  va_end(ap);
}


void
debug_printf(const int debug_level, const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  vfprintf(my_dbg, fmt, ap);
  va_end(ap);
}

void
error_printf(const char *fmt, ...)
{
  va_list ap;

  // Output to error stream.
  va_start(ap, fmt);
  vfprintf(my_err, fmt, ap);
  va_end(ap);

  // Also output to debug stream.
  va_start(ap, fmt);
  vfprintf(my_dbg, fmt, ap);
  va_end(ap);
}

/**
 * @brief Output resource usage statistics
 */
void
output_procinfo(const output_type ptype, const int level)
{
  int status;
  struct rusage my_rusage;
  int uhours;
  int uminutes;
  int useconds;
  int shours;
  int sminutes;
  int sseconds;
  FILE *pipe;
  char ps_command[1024];
  pid_t pid;
  unsigned long rss_mem_usage;
  unsigned long vsz_mem_usage;
  unsigned long total_mem_usage;

  // Find out how much CPU time process has used.
  status = getrusage(RUSAGE_SELF, &my_rusage);
  assert(status != -1);

  // Calculate hours:minutes:seconds for process 'user' time.
  useconds = my_rusage.ru_utime.tv_sec;
  uhours = useconds / 3600;
  useconds -= uhours * 3600;
  uminutes = useconds / 60;
  useconds -= uminutes * 60;

  // Calculate hours:minutes:seconds for process 'system' time.
  sseconds = my_rusage.ru_stime.tv_sec;
  shours = sseconds / 3600;
  sseconds -= shours * 3600;
  sminutes = sseconds / 60;
  sseconds -= sminutes * 60;

  // Note process ID.
  pid = getpid();

  // Run ps command to get resident memory usage.
  sprintf(ps_command, "ps -o rss -p %d 2>/dev/null | grep -v RSS", pid);
  pipe = popen(ps_command, "r");
  if (NULL == pipe)
  {
    debug_printf(0, "Error running ps command to get memory usage (RSS)\n");
    rss_mem_usage = 0;
  }
  else
  {
    char buf[1024];
    char *ptr;

    fgets(buf, sizeof(buf), pipe);

    // Skip leading spaces
    ptr = buf;
    while (('\0' != *ptr) && (' ' == *ptr))
    {
      ptr += 1;
    }

    rss_mem_usage = strtoul(ptr, NULL, 10);
    pclose(pipe);
  }

  // Run ps command to get virtual memory usage.
  sprintf(ps_command, "ps -o vsz -p %d | grep -v VSZ", pid);
  pipe = popen(ps_command, "r");
  if (NULL == pipe)
  {
    debug_printf(0, "Error running ps command to get memory usage (VSZ)\n");
    vsz_mem_usage = 0;
  }
  else
  {
    char buf[1024];
    char *ptr;

    fgets(buf, sizeof(buf), pipe);

    // Skip leading spaces
    ptr = buf;
    while (('\0' != *ptr) && (' ' == *ptr))
    {
      ptr += 1;
    }

    vsz_mem_usage = strtoul(ptr, NULL, 10);
    pclose(pipe);
  }


  // Round memory usage to nearest megabyte.
  total_mem_usage = (rss_mem_usage + vsz_mem_usage + 512) / 1024;
  rss_mem_usage = (rss_mem_usage + 512) / 1024;
  vsz_mem_usage = (vsz_mem_usage + 512) / 1024;

  custom_printf(ptype, level,
                "usertime:%02d:%02d:%02d.%d systime:%02d:%02d:%02d.%d "
                "RSS:%luMB; VSZ:%luMB; Total memory:%luMB\n",
                uhours, uminutes, useconds,
                (int)(my_rusage.ru_utime.tv_usec / 1000),
                shours, sminutes, sseconds,
                (int)(my_rusage.ru_utime.tv_usec / 1000),
                rss_mem_usage, vsz_mem_usage, total_mem_usage);
}

