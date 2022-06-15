#ifndef UTIL_LOG_H_
#define UTIL_LOG_H_

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <memory>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
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
  LogMessage() : LogMessage(LogLevel::kInfo, "") {}
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
  virtual void Write(const char* str) = 0;
  void WriteLine(const LogMessage& msg);

 protected:
  LogWriter();
  LogFormatterPtrU fmt_;
};


class LogConsoleWriter : public LogWriter {
 public:
  void Write(const LogMessage& msg) override;
  void Write(const char* str) override;
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
  void Emit(bool single_line, LogLevel lv, const char* fmt, Args&&... args) {
    LogMessage msg{};
    if constexpr (sizeof...(args) > 0) {
      std::snprintf(buf_, kBufLen, fmt, std::forward<Args>(args)...);
      msg = LogMessage{ lv, buf_ };
    } else {
      msg = LogMessage(lv, fmt);
    }
    for (const auto& [f, w] : writers_) {
      if (f->Check(msg)) {
        single_line ? w->WriteLine(msg) : w->Write(msg.msg_.c_str());
      }
    }
  }

  static constexpr size_t kBufLen = 1024;

 private:
  Logger();

  char buf_[kBufLen];
  std::vector<std::tuple<LogFilterPtrU, LogWriterPtrU>> writers_;
};


#ifdef FOR_TEST
#define LOG_DEBUG(...) halo_pm::Logger::GetInstance().Emit(true, halo_pm::LogLevel::kDebug, __VA_ARGS__)
#else
#define LOG_DEBUG(...)
#endif
#define LOG_VERBOSE(...) halo_pm::Logger::GetInstance().Emit(true, halo_pm::LogLevel::kVerbose, __VA_ARGS__)
#define LOG_INFO(...) halo_pm::Logger::GetInstance().Emit(true, halo_pm::LogLevel::kInfo, __VA_ARGS__)
#define LOG_WARNING(...) halo_pm::Logger::GetInstance().Emit(true, halo_pm::LogLevel::kWarning, __VA_ARGS__)
#define LOG_ERROR(...) halo_pm::Logger::GetInstance().Emit(true, halo_pm::LogLevel::kError, __VA_ARGS__)

#ifdef FOR_TEST
#define LOG_DEBUG_C(...) halo_pm::Logger::GetInstance().Emit(false, halo_pm::LogLevel::kDebug, __VA_ARGS__)
#else
#define LOG_DEBUG_C(...)
#endif
#define LOG_VERBOSE_C(...) halo_pm::Logger::GetInstance().Emit(false, halo_pm::LogLevel::kVerbose, __VA_ARGS__)
#define LOG_INFO_C(...) halo_pm::Logger::GetInstance().Emit(false, halo_pm::LogLevel::kInfo, __VA_ARGS__)
#define LOG_WARNING_C(...) halo_pm::Logger::GetInstance().Emit(false, halo_pm::LogLevel::kWarning, __VA_ARGS__)
#define LOG_ERROR_C(...) halo_pm::Logger::GetInstance().Emit(false, halo_pm::LogLevel::kError, __VA_ARGS__)

}  // namespace halo_pm

#endif  // UTIL_LOG_H_
