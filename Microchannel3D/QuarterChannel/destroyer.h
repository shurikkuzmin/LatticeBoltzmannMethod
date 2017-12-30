#ifndef DESTROYER_H
#define DESTROYER_H

#include <list>

template<typename T>
class Destroyer
{
	typedef std::list<T*> List;
	typedef typename List::iterator ListIterator;

public:

	Destroyer()
	{
	}

	~Destroyer()
	{
		for(ListIterator k = objects.begin(); k != objects.end(); ++k)
		{
			delete *k;
		}
	}

	T * take_ownership(T * object)
	{
		objects.push_back(object);
		return object;
	}

private:

	Destroyer(const Destroyer &);
	Destroyer & operator=(const Destroyer &);

	List objects;
};

#endif
