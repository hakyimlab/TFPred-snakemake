


library(inline)

wc.code <- "
uintmax_t linect = 0; 
uintmax_t tlinect = 0;

int fd, len;
u_char *p;

struct statfs fsb;

static off_t buf_size = SMALL_BUF_SIZE;
static u_char small_buf[SMALL_BUF_SIZE];
static u_char *buf = small_buf;

PROTECT(f = AS_CHARACTER(f));

if ((fd = open(CHAR(STRING_ELT(f, 0)), O_RDONLY, 0)) >= 0) {

  if (fstatfs(fd, &fsb)) {
    fsb.f_iosize = SMALL_BUF_SIZE;
  }

  if (fsb.f_iosize != buf_size) {
    if (buf != small_buf) {
      free(buf);
    }
    if (fsb.f_iosize == SMALL_BUF_SIZE || !(buf = malloc(fsb.f_iosize))) {
      buf = small_buf;
      buf_size = SMALL_BUF_SIZE;
    } else {
      buf_size = fsb.f_iosize;
    }
  }

  while ((len = read(fd, buf, buf_size))) {

    if (len == -1) {
      (void)close(fd);
      break;
    }

    for (p = buf; len--; ++p)
      if (*p == '\\n')
        ++linect;
  }

  tlinect += linect;

  (void)close(fd);

}
SEXP result;
PROTECT(result = NEW_INTEGER(1));
INTEGER(result)[0] = tlinect;
UNPROTECT(2);
return(result);
";

setCMethod("wc",
           signature(f="character"), 
           wc.code,
           includes=c("#include <stdlib.h>", 
                      "#include <stdio.h>",
                      "#include <sys/param.h>",
                      "#include <sys/mount.h>",
                      "#include <sys/stat.h>",
                      "#include <ctype.h>",
                      "#include <err.h>",
                      "#include <errno.h>",
                      "#include <fcntl.h>",
                      "#include <locale.h>",
                      "#include <stdint.h>",
                      "#include <string.h>",
                      "#include <unistd.h>",
                      "#include <wchar.h>",
                      "#include <wctype.h>",
                      "#define SMALL_BUF_SIZE (1024 * 8)"),
           language="C",
           convention=".Call")

wc("FULLPATHTOFILE")