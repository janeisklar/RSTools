#include "unix.hpp"

namespace rstools {
namespace batch {
namespace execution {
    
void Unix::_parseParams(int argc, char * argv[])
{
    this->executionSuccessful = true;
}
    
void Unix::_init() {}

void Unix::_run()
{
    if ( system(this->getTask()->getCmd()) != 0 ) {
        this->executionSuccessful = false;
    }
}

void Unix::destroy() {}

bool Unix::isEverythingFine()
{
    return this->executionSuccessful;
}

}}} // namespace rstools::batch::execution
