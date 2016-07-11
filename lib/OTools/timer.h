/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#ifndef _TIMER_H
#define _TIMER_H

#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <string>

class timer {
protected:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_timing_clock, prev_timing_clock;

public:
	timer () {
		start_timing_clock = std::chrono::high_resolution_clock::now();
	}

	~timer() {
	}

	void clock() {
		prev_timing_clock = std::chrono::high_resolution_clock::now();
	}

	unsigned int rel_time() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - prev_timing_clock).count();
	}

	unsigned int abs_time() {
		return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_timing_clock).count();
	}

	std::string date() {
		char buffer[256];
		std::time_t t = std::time(NULL);
	    std::tm tm = *std::localtime(&t);
	    std::strftime(buffer, 256, "%d/%m/%Y - %X", &tm);
	    return string(buffer);
	}

	/* More elegant GCC5 version
	std::string date() {
		auto now = std::chrono::system_clock::now();
		auto in_time_t = std::chrono::system_clock::to_time_t(now);
		std::stringstream ss;
	    ss << std::put_time(std::localtime(&in_time_t), "%d/%m/%Y - %X");
	    return ss.str();
	}
	*/

};

#endif
