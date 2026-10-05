#ifndef _CONFIGURATION_OPTIONS_HPP_
#define _CONFIGURATION_OPTIONS_HPP_

#include <string>
#include <map>

enum class VerboseOption {SILENT, DISABLE_WARNING, VERBOSE};
enum class ModelOption {HORN_SCHUNCK, CORNELIUS_KANADE};

inline const std::map<std::string, VerboseOption> mapper_verbose_option {
    {"silent", VerboseOption::SILENT},
    {"disable-warnings", VerboseOption::DISABLE_WARNING},
    {"verbose", VerboseOption::VERBOSE}
};

inline const std::map<std::string, ModelOption> mapper_model_option {
    {"horn-schunck", ModelOption::HORN_SCHUNCK},
    {"cornelius-kanade", ModelOption::CORNELIUS_KANADE}
};

#endif