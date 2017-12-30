#ifndef __MIXER3D_SINGLETON__
#define __MIXER3D_SINGLETON__

template<class T> class Singleton
{
	static T * instance;
	static T * create_instance();

public:

	Singleton()
	{
		instance = create_instance();
	}

	static T & get_instance()
	{
		return *instance;
	}

	virtual ~Singleton()
	{
		delete instance;
		instance = 0;
	}
};

template<class T> class LazySingleton
{
	static T * instance;
	static T * create_instance();

public:

	static T & get_instance()
	{
		if(instance == 0)
			instance = create_instance();
		return *instance;
	}

	virtual ~LazySingleton()
	{
		delete instance;
		instance = 0;
	}
};

#endif
