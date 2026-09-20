port <- as.integer(Sys.getenv("COTRA_PORT", "3838"))
con <- try(
  socketConnection(
    host = "127.0.0.1",
    port = port,
    open = "r+b",
    timeout = 3
  ),
  silent = TRUE
)

if (inherits(con, "try-error")) {
  quit(status = 1)
}

close(con)
quit(status = 0)
