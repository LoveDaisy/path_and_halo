#ifndef UTIL_LOG_H_
#define UTIL_LOG_H_

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <memory>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace halo_pm {

class LogFilter;
class LogFormatter;
class LogWriter;

using LogFilterPtrU = std::unique_ptr<LogFilter>;
using LogFormatterPtrU = std::unique_ptr<LogFormatter>;
using LogWriterPtrU = std::unique_ptr<LogWriter>;

enum class LogLevel : uint8_t {
  kDebug,
  kVerbose,
  kInfo,
  kWarning,
  kError,
};

inline bool operator==(LogLevel l, LogLevel r) {
  return static_cast<uint8_t>(l) == static_cast<uint8_t>(r);
}

inline bool operator<(LogLevel l, LogLevel r) {
  return static_cast<uint8_t>(l) < static_cast<uint8_t>(r);
}


struct LogMessage {
  LogLevel lv_;
  std::chrono::time_point<std::chrono::system_clock> t_;
  std::thread::id thread_id_;
  std::string msg_;

  LogMessage(LogLevel lv, const char* msg);
};


// LogFormatter
class LogFormatter {
 public:
  LogFormatter() = default;
  LogFormatter(const LogFormatter&) = default;
  LogFormatter(LogFormatter&&) = default;
  virtual ~LogFormatter() = default;

  LogFormatter& operator=(const LogFormatter&) = default;
  LogFormatter& operator=(LogFormatter&&) = default;

  virtual const char* Format(const LogMessage& msg) = 0;

  static constexpr size_t kBufLen = 1024 * 32;

 protected:
  char buf_[kBufLen];
};


class SimpleLogFormatter : public LogFormatter {
 public:
  const char* Format(const LogMessage& msg) override;
};


// LogWriter
class LogWriter {
 public:
  LogWriter(const LogWriter&) = delete;
  LogWriter(LogWriter&&) = default;
  virtual ~LogWriter() = default;

  LogWriter& operator=(const LogWriter&) = delete;
  LogWriter& operator=(LogWriter&&) = default;

  void SetFormatter(LogFormatterPtrU fmt);
  virtual void Write(const LogMessage& msg) = 0;

 protected:
  LogWriter();
  LogFormatterPtrU fmt_;
};


class LogConsoleWriter : public LogWriter {
 public:
  void Write(const LogMessage& msg) override;
};


// LogFilter
class LogFilter {
 public:
  LogFilter() = default;
  LogFilter(const LogFilter&) = default;
  LogFilter(LogFilter&&) = default;
  virtual ~LogFilter() = default;

  LogFilter& operator=(const LogFilter&) = default;
  LogFilter& operator=(LogFilter&&) = default;

  virtual bool Check(const LogMessage& msg) = 0;
};


class MinLogLevelFilter : public LogFilter {
 public:
  MinLogLevelFilter(LogLevel min_lv);
  bool Check(const LogMessage& msg) override;

 private:
  LogLevel min_lv_;
};


// Logger
class Logger {
 public:
  static Logger& GetInstance();

  template <class... Args>
  void Emit(LogLevel lv, const char* fmt, Args&&... args) {
    std::snprintf(buf_, kBufLen, fmt, std::forward<Args>(args)...);
    LogMessage msg(lv, buf_);
    for (const auto& [f, w] : writers_) {
      if (f->Check(msg)) {
        w->Write(msg);
      }
    }
  }

  static constexpr size_t kBufLen = 1024;

 private:
  Logger();

  char buf_[kBufLen];
  std::vector<std::tuple<LogFilterPtrU, LogWriterPtrU>> writers_;
};


#define LOG_DEBUG(...) halo_pm::Logger::GetInstance().Emit(halo_pm::LogLevel::kDebug, __VA_ARGS__)
#define LOG_VERBOSE(...) halo_pm::Logger::GetInstance().Emit(halo_pm::LogLevel::kVerbose, __VA_ARGS__)
#define LOG_INFO(...) halo_pm::Logger::GetInstance().Emit(halo_pm::LogLevel::kInfo, __VA_ARGS__)
#define LOG_WARNING(...) halo_pm::Logger::GetInstance().Emit(halo_pm::LogLevel::kWarning, __VA_ARGS__)
#define LOG_ERROR(...) halo_pm::Logger::GetInstance().Emit(halo_pm::LogLevel::kError, __VA_ARGS__)

}  // namespace halo_pm

#endif  // UTIL_LOG_H_
