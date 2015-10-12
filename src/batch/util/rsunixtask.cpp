#include "rsunixtask.hpp"

namespace rstools {
namespace batch {
namespace util {

RSUnixTask::RSUnixTask(const char* code, const char* name) : RSTask(code, name)
{}

bool RSUnixTask::hasInputNiftiHeaderInformation()
{
    return this->inputNiftiHeaderInformation != NULL;
}

rsNiftiExtendedHeaderInformation* RSUnixTask::getInputNiftiHeaderInformation()
{
    return this->inputNiftiHeaderInformation;
}

void RSUnixTask::setInputNiftiHeaderInformation(rsNiftiExtendedHeaderInformation* info)
{
    this->inputNiftiHeaderInformation = info;
}

char * RSUnixTask::getTempDirectoryPath()
{
    return this->tempDirectoryPath;
}

void RSUnixTask::setTempDirectoryPath(char *path)
{
    this->tempDirectoryPath = path;
}

}}} // namespace rstools::batch::util
