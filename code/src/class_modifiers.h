#ifndef class_modifiers_h
#define class_modifiers_h

class uncopyable {
public:
	uncopyable() {}
private:
	uncopyable(const uncopyable&);
	uncopyable& operator=(const uncopyable&);	
};

#endif //class_modifiers_h
