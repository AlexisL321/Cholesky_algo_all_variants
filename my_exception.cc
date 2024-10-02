#include "my_exception.h"

using namespace std;

//constructor
My_exception::My_exception(const string& message) : m_message(message){}

//override what() to return custom message
const char* My_exception::what() const noexcept {
	return m_message.c_str();
}
