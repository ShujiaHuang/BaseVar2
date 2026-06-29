/**
 * @file dump.h
 *
 * @brief Dump tool for inspecting binary batchfile (.bbf) and binary index (.bbi) files.
 *
 * Usage:
 *   basevar dump <file.bbf>  [--header] [--region chr:start-end] [--limit N] [-v]
 *   basevar dump <file.bbi>  [--entries]
 *
 * @author Shujia Huang
 * @date 2026
 */
#ifndef __INCLUDE_BASEVAR_DUMP_H__
#define __INCLUDE_BASEVAR_DUMP_H__

namespace basevar {
namespace dump {

/**
 * @brief Dump subcommand entry point.
 *
 *  argv layout: argv[0] == "dump", then options + file path.
 *
 * @param argc  Argument count.
 * @param argv  Argument vector.
 * @return Exit code (0 on success, non-zero on error).
 */
int dump_runner(int argc, char* argv[]);

}  // namespace dump
}  // namespace basevar

#endif
