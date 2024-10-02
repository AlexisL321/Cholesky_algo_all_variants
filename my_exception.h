#ifndef MY_EXCEPTION_H
#define MY_EXCEPTION_H

#include <stdexcept>
#include <string>

using namespace std;

class My_exception : public exception{
	public:
		My_exception(const string& message);

		const char* what() const noexcept override;

	private:
		string m_message;
};

#endif
