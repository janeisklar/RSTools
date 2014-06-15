#ifndef rstools_rstask_hpp
#define rstools_rstask_hpp

#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <locale>
#include <vector>

using namespace std;

typedef struct {
	char *key;
	char *value;
} rsArgument;

class RSTask {
	
	public:
		static const short TASK_UNKNOWN           = -1;
		static const short TASK_UNIX              = 0;
		static const short TASK_RSTIMECOURSE      = 1;
		static const short TASK_RSREGRESSION      = 2;
		static const short TASK_RSBANDPASS        = 3;
		static const short TASK_RSMOTIONSCRUBBING = 4;
		static const short TASK_RSCORRELATION     = 5;
		static const short TASK_RSROI             = 6;
		
		RSTask(short);
		~RSTask();
		
		short getTask();
		void setDescription(char*);
		char* getDescription();
		void addArgument(rsArgument*);
		vector<rsArgument*> getArguments();
		char* getOutputPath();
		void setOutputPath(char*);
		void setShowOutput(bool showOutput);
		bool shouldShowOutput();
		char** getCallString(int *argc);
		char const * getName();
		void setCmd(char* cmd);
		char* getCmd();
		
		static short getTaskFromName(char*);
		static char const* getNameForTask(short);
	
	protected:
		short task;
		char *description;
		vector<rsArgument*> arguments;
		char *outputPath;
		char *cmd;
		bool showOutput;
		
		static int strcicmp(char const *, char const *);
};

#endif
