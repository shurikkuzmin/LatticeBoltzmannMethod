#ifndef __PARAMETERS_LIST_H__
#define __PARAMETERS_LIST_H__

#include <vector>
#include <map>
#include <string>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

class ParamsList
{
public:

	typedef std::vector<std::string>      Values;
	typedef std::map<std::string, Values> Params;
	typedef Params::iterator              Iter;
	typedef Params::const_iterator        ConstIter;

	class Parameter
	{
	public:

		explicit Parameter(ConstIter iter) : iter(iter) {}

		size_t size()const
		{
			return iter->second.size();
		}

		template<typename T>
		T value(size_t index = 0)const
		{
			return boost::lexical_cast<T>(this->value(index));
		}

// 		template<typename T>
// 		T& value(size_t index = 0)
// 		{
// 			return boost::lexical_cast<T>(this->value(index));
// 		}

		std::string value(size_t index = 0)const
		{
			if(index >= this->size())
			{
				throw std::runtime_error("Parameter value index is out of bounds.");
			}

			return iter->second[index];
		}

	private:
		ConstIter iter;
	};

	ParamsList()
	{}

	ParamsList(const ParamsList & p)
		:	params(p.params)
	{}

	ParamsList & operator=(const ParamsList & p)
	{
		if(&p != this)
		{
			this->params = p.params;
		}
		return *this;
	}

	bool has(const std::string & name)const
	{
		ConstIter k = params.find(name);
		return k != params.end();
	}

	const Parameter operator()(const std::string & name)const
	{
		ConstIter k = params.find(name);
		if(k == params.end())
		{
			throw std::runtime_error("In ParameterList::get(): parameter does not exist: " + name);
		}

		return Parameter(k);
	}

	
	template<typename T>
	void change(const std::string & name, const T& value)
	{
		if(has(name))
			params[name].front()=boost::lexical_cast<std::string>(value);
	}


	void add(const std::string & name, const Values & values)
	{
		std::pair<Iter, bool> result = params.insert(std::make_pair(name, values));
		if(!result.second)
		{
			throw std::runtime_error(std::string("Parameter with this name already exists in ParameterList: ") + name);
		}
	}

	void add(const std::string & name, const std::string & value)
	{
		Values values;
		values.push_back(value);
		this->add(name, values);
	}

	template<typename T>
	void add(const std::string & name, const T & value)
	{
		Values values;
		values.push_back(boost::lexical_cast<std::string>(value));
		this->add(name, values);
	}

private:

	Params params;
};

#endif
