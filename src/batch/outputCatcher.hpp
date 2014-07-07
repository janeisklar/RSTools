#ifndef rstools_batcb_outputcatcher_hpp
#define rstools_batcb_outputcatcher_hpp

#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <cstring>

namespace rstools {
namespace batch {

class OutputCatcher
{
public:

	void initialize()
	{
		if (m_init) {
			return;
		}
		m_pipe[READ] = 0;
        m_pipe[WRITE] = 0;
        if (pipe(m_pipe) == -1) {
            return;
		}
        m_oldStdOut = dup(fileno(stdout));
        m_oldStdErr = dup(fileno(stderr));
        if (m_oldStdOut == -1 || m_oldStdErr == -1) {
            return;
		}
		
		oldStdout = fdopen(m_oldStdOut, "a");

        m_init = true;
	}
	
	void tearDown()
	{
		if ( ! m_init ) {
			return;
		}
        if (m_oldStdOut > 0) {
            close(m_oldStdOut);
		}
        if (m_oldStdErr > 0) {
            close(m_oldStdErr);
		}
        if (m_pipe[READ] > 0) {
            close(m_pipe[READ]);
		}
        if (m_pipe[WRITE] > 0) {
            close(m_pipe[WRITE]);
		}
		
		if ( oldStdout != NULL ) {
			fclose(oldStdout);
		}
		
		m_init = false;
	}

    OutputCatcher()
    {
		m_capturing = false;
		m_init = false;
		m_oldStdOut = 0;
		m_oldStdErr = 0;
    }

    ~OutputCatcher()
    {
		tearDown();
    }

    void beginCapture()
    {
        if (!m_init)
            initialize();
        if (m_capturing)
            endCapture();
        fflush(stdout);
        fflush(stderr);
        dup2(m_pipe[WRITE], fileno(stdout));
        dup2(m_pipe[WRITE], fileno(stderr));
        m_capturing = true;
    }

    bool endCapture()
    {
        if (!m_init) {
            return false;
		}
        if (!m_capturing) {
			tearDown();
            return false;
		}
		fprintf(stdout, "\n"); // print \n so that our pipe will receive something even if nothing is being sent
        fflush(stdout);
        fflush(stderr);
        dup2(m_oldStdOut, fileno(stdout));
        dup2(m_oldStdErr, fileno(stderr));
        m_captured.clear();

        char* buf;
        const int bufSize = 32768;
        buf = (char*)malloc((size_t)bufSize*sizeof(char));
        int bytesRead = 0;

        do {
            bytesRead = read(m_pipe[READ], buf, bufSize-1);
			buf[bytesRead] = '\0';
            m_captured += buf;
		} while(bytesRead == bufSize);

		tearDown();
		
        return true;
    }

    std::string getCapture() const
    {
		// ignore first character, because it just ensures that read() in endCapture() actually receives some data
        return m_captured.substr(0, m_captured.length()-1);
    }

	FILE* getStdout()
	{
		return oldStdout;
	}

private:
    enum PIPES { READ, WRITE };
    int m_pipe[2];
    int m_oldStdOut;
    int m_oldStdErr;
    bool m_capturing;
    bool m_init;
    std::string m_captured;
	FILE *oldStdout;
};

}}

#endif
