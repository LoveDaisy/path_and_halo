#include "util/log.hpp"

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <memory>
#include <ratio>
#include <thread>
#include <utility>

namespace halo_pm {

// LogMessage
LogMessage::LogMessage(LogLevel lv, const char* msg)
    : lv_(lv), t_(std::chrono::system_clock::now()), thread_id_(std::this_thread::get_id()), msg_(msg) {}


constexpr LogLevel kDefaultFilterLevel =
#ifdef FOR_TEST
    LogLevel::kDebug;
#else
    LogLevel::kInfo;
#endif


// Logger
Logger::Logger() : buf_{} {
  writers_.emplace_back(LogFilterPtrU{ new MinLogLevelFilter(kDefaultFilterLevel) },
                        LogWriterPtrU{ new LogConsoleWriter });
};

Logger& Logger::GetInstance() {
  static Logger logger;
  return logger;
}


// LogFormatter
const char* SimpleLogFormatter::Format(const LogMessage& msg) {
  size_t thread_id = std::hash<std::thread::id>{}(msg.thread_id_);

  const char* lv_str = "NONE";
  switch (msg.lv_) {
    case LogLevel::kDebug:
      lv_str = "DEBUG";
      break;
    case LogLevel::kVerbose:
      lv_str = "VERBOSE";
      break;
    case LogLevel::kInfo:
      lv_str = "INFO";
      break;
    case LogLevel::kWarning:
      lv_str = "WARNING";
      break;
    case LogLevel::kError:
      lv_str = "ERROR";
      break;
  }

  size_t offset = 0;
  offset += std::snprintf(buf_, kBufLen, "%06zu[%s]", thread_id % 1000000, lv_str);

  auto t = std::chrono::system_clock::to_time_t(msg.t_);
  offset += std::strftime(buf_ + offset, kBufLen, "%H:%M:%S.", std::localtime(&t));

  auto t0 = std::chrono::floor<std::chrono::seconds>(msg.t_);
  auto t_ms = std::chrono::duration<double, std::milli>(msg.t_ - t0);
  std::snprintf(buf_ + offset, kBufLen, "%03d %s", static_cast<int>(t_ms.count()), msg.msg_.c_str());
  return buf_;
}


// LogWriter
LogWriter::LogWriter() : fmt_(new SimpleLogFormatter) {}

void LogWriter::SetFormatter(LogFormatterPtrU fmt) {
  fmt_ = std::move(fmt);
}

void LogConsoleWriter::Write(const LogMessage& msg) {
  const char* msg_str = msg.msg_.c_str();
  if (fmt_) {
    msg_str = fmt_->Format(msg);
  }
  if (msg.lv_ < LogLevel::kWarning) {
    std::printf("%s\n", msg_str);
  } else {
    std::fprintf(stderr, "%s\n", msg_str);
  }
}


// LogFilter
MinLogLevelFilter::MinLogLevelFilter(LogLevel min_lv) : min_lv_(min_lv) {}

bool MinLogLevelFilter::Check(const LogMessage& msg) {
  return msg.lv_ >= min_lv_;
}

}  // namespace halo_pm
