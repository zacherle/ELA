
//---------------------------------------------------------------------
// Supporting datetime and timedelta classes (equivalent to datetime_module)
//---------------------------------------------------------------------
class timedelta; // forward declaration

class datetime {
private:
    std::chrono::system_clock::time_point tp;
public:
    // Default constructor
    datetime() {
        tp = std::chrono::system_clock::now();
    }
    // Constructor from year only (e.g., datetime(3000))
    explicit datetime(int year) {
        std::tm t = {};
        t.tm_year = year - 1900;
        t.tm_mon = 0;
        t.tm_mday = 1;
        t.tm_hour = 0;
        t.tm_min = 0;
        t.tm_sec = 0;
        std::time_t time_tt = std::mktime(&t);
        tp = std::chrono::system_clock::from_time_t(time_tt);
    }
    // Constructor from full date and time (seconds and millisecond)
    datetime(int year, int month, int day, int hour, int minute, int second, int msec=0) {
        std::tm t = {};
        t.tm_year = year - 1900;
        t.tm_mon = month - 1;
        t.tm_mday = day;
        t.tm_hour = hour;
        t.tm_min = minute;
        t.tm_sec = second;
        std::time_t time_tt = std::mktime(&t);
        tp = std::chrono::system_clock::from_time_t(time_tt) + std::chrono::milliseconds(msec);
    }
    // Copy constructor
    datetime(const datetime &other) {
        tp = other.tp;
    }
    // Overloaded assignment operator
    datetime& operator=(const datetime &other) {
        if(this != &other) {
            tp = other.tp;
        }
        return *this;
    }
    // Addition of timedelta
    datetime operator+(const timedelta &td) const;
    // Subtraction of two datetime objects gives a timedelta
    timedelta operator-(const datetime &other) const;
    // Returns ISO formatted date string
    std::string isoformat() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()) % 1000;
        std::ostringstream oss;
        oss << std::put_time(ptm, "%Y-%m-%dT%H:%M:%S") << "." << std::setw(3) << std::setfill('0') << ms.count();
        return oss.str();
    }
    // Getters for datetime components
    int getYear() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_year + 1900;
    }
    int getMonth() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_mon + 1;
    }
    int getDay() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_mday;
    }
    int getHour() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_hour;
    }
    int getMinute() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_min;
    }
    int getSecond() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_sec;
    }
    int getMillisecond() const {
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()) % 1000;
        return static_cast<int>(ms.count());
    }
    // Comparison operator (less than)
    bool operator<(const datetime &other) const {
        return tp < other.tp;
    }
};

class timedelta {
private:
    // Store duration in milliseconds
    std::chrono::milliseconds duration;
public:
    // Constructor from milliseconds value via parameters (e.g., timedelta(0,0,0,0, rtime_m))
    timedelta(int days, int hours, int minutes, int seconds, int millis) {
        long long total_ms = millis;
        total_ms += static_cast<long long>(seconds) * 1000;
        total_ms += static_cast<long long>(minutes) * 60 * 1000;
        total_ms += static_cast<long long>(hours) * 3600 * 1000;
        total_ms += static_cast<long long>(days) * 24 * 3600 * 1000;
        duration = std::chrono::milliseconds(total_ms);
    }
    // Default constructor
    timedelta() : duration(0) {}
    // Getter for total seconds as double
    double total_seconds() const {
        return duration.count() / 1000.0;
    }
    // Allow addition of two timedeltas
    timedelta operator+(const timedelta &other) const {
        timedelta temp;
        temp.duration = duration + other.duration;
        return temp;
    }
    // Allow subtraction of two timedeltas
    timedelta operator-(const timedelta &other) const {
        timedelta temp;
        temp.duration = duration - other.duration;
        return temp;
    }
    // Expose the internal duration for use in datetime arithmetic
    const std::chrono::milliseconds& getDuration() const {
        return duration;
    }
};

// Definition of datetime + timedelta operator
datetime datetime::operator+(const timedelta &td) const {
    datetime temp(*this);
    temp = datetime(); // reset temp to force assignment below
    temp.tp = this->tp + td.getDuration();
    return temp;
}

// Definition of datetime - datetime operator to yield timedelta
timedelta datetime::operator-(const datetime &other) const {
    timedelta td;
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(this->tp - other.tp);
    // Directly assign the computed difference in milliseconds
    td = timedelta(0,0,0,0, static_cast<int>(diff.count()));
    return td;
}

