// 'MAGIC' parallelization module written by Nils Schöneberg, adapted from the one created by Thomas Tram

//To be used for declaring which arguments a parallel clause 'captures' from the outside
//This would be equivalent to the 'firstprivate' in OPENMP
#define with_arguments(...) __VA_ARGS__
//To declare a list of arguments a parallel clause only has internally
//(necessary since comma-separated lists are interpreted as separate arguments
// except when surrounded by backets as forced by this macro)
#define declare_list_of_variables_inside_parallel_region(...) __VA_ARGS__

// To be called WITHIN a given parallel loop. First argument is the parameters to 'capture',
// corresponding to the 'firstprivate' in OPENMP.
// Usually either '=' to capture all necessary variables, or a list of arguments
// declared within the 'with_arguments' macro for a more precise list
// The second argument is the actual code to execute. That code should be formulated like
// the body of a normal function, e.g. return _SUCCESS_ at the end.
// 'private' arguments in OPENMP can be emulated by just making them local to the function.
// Special care: A list of private arguments needs to be surrounded by the special macro
// 'declare_list_of_variables-inside_parallel_region' due to the peculiarities
// of the C++ preprocessor macros. See examples e.g.
// in source/perturbations.c, source/lensing.c, or tools/hypershperical.c
#define class_run_parallel(arg1, arg2) future_output.push_back(task_system.AsyncTask([arg1] () {arg2}));
// The mutable version allows one to change variables outside the scope of the parallel region
// Be careful, this is very dangerous, and you should be sure that you don't access any
// variables you edit during the parallel region outside of it, as that value will be random
#define class_run_parallel_mutable(arg1, arg2) future_output.push_back(task_system.AsyncTask([arg1] () mutable {arg2}));

// To be called ONLY ONCE without arguments before the intended parallel loop(s).
#define class_setup_parallel()                    \
Tools::TaskSystem task_system;                    \
std::vector<std::future<int>> future_output;

// To be called ONLY ONCE without arguments before the intended parallel loop(s).
#define class_setup_parallel_optional(is_multi_threaded) \
Tools::TaskSystem task_system{ (is_multi_threaded) ? Tools::TaskSystem::GetNumThreads() : 1 }; \
std::vector<std::future<int>> future_output;

// To be called without arguments AFTER ANY parallel loop in order to actually execute the jobs.
// NEEDS TO BE CALLED BEFORE USING THE RESULTS!
#define class_finish_parallel()                   \
for (std::future<int>& future : future_output) {  \
  if(future.get()!=_SUCCESS_) return _FAILURE_;   \
}                                                 \
future_output.clear();

//
//  thread_pool.h
//  ppCLASS
//
//  Created by Thomas Tram on 02/03/2020.
//  Copyright 2020 Aarhus University. All rights reserved.
//
#ifndef THREAD_POOL_H
#define THREAD_POOL_H
#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace Tools {

class NotificationQueue {
public:
  bool TryPop(std::function<void()>& x) {
    std::unique_lock<std::mutex> lock(mutex_, std::try_to_lock);
    if (!lock || queue_.empty()) {
      return false;
    }
    x = std::move(queue_.front());
    queue_.pop_front();
    return true;
  }

  bool Pop(std::function<void()>& x) {
    std::unique_lock<std::mutex> lock(mutex_);
    while (queue_.empty() && !done_) {
      ready_.wait(lock);
    }
    if (queue_.empty()) {
      return false;
    }
    x = std::move(queue_.front());
    queue_.pop_front();
    return true;
  }

  template<typename F>
  bool TryPush(F&& f) {
    {
      std::unique_lock<std::mutex> lock(mutex_, std::try_to_lock);
      if (!lock) {
        return false;
      }
      queue_.emplace_back(std::forward<F>(f));
    }
    ready_.notify_one();
    return true;
  }

  template<typename F>
  void Push(F&& f) {
    {
      std::unique_lock<std::mutex> lock(mutex_);
      queue_.emplace_back(std::forward<F>(f));
    }
    ready_.notify_one();
  }

  void Done() {
    {
      std::unique_lock<std::mutex> lock(mutex_);
      done_ = true;
    }
    ready_.notify_all();
  }
private:
  std::deque<std::function<void()>> queue_;
  bool done_ = false;
  std::mutex mutex_;
  std::condition_variable ready_;
};

class TaskSystem {
public:
  TaskSystem(unsigned int count = GetNumThreads())
  : count_(count)
  , index_(0)
  , queues_{count_} {
    for (unsigned int n = 0; n < count_; ++n) {
      threads_.emplace_back([&, n]{ Run(n); });
    }
  }

  ~TaskSystem() {
    for (auto& e : queues_) e.Done();
    for (auto& e : threads_) e.join();
  }

  static unsigned int GetNumThreads() {
    unsigned int number_of_threads = std::thread::hardware_concurrency();
    for (const std::string& env_var_name : {"OMP_NUM_THREADS", "SLURM_CPUS_PER_TASK"}) {
      if (char* s = std::getenv(env_var_name.c_str())) {
        int threads = std::atoi(s);
        if ((threads > 0) && (threads <= 8192)) {
          number_of_threads = threads;
          break;
        }
      }
    }
    return number_of_threads;
  }

  template <typename F>
  void Async(F&& f) {
    auto i = index_++;
    for (unsigned int n = 0; n < count_; ++n) {
      if (queues_[(i + n) % count_].TryPush(std::forward<F>(f))) {
        return;
      }
    }
    queues_[i % count_].Push(std::forward<F>(f));
  }

  template<typename F>
  std::future<typename std::result_of<F()>::type> AsyncTask(F&& f) {
    using return_type = typename std::result_of<F()>::type;
    auto task = std::make_shared<std::packaged_task<return_type()>>(f);
    std::future<return_type> res = task->get_future();

    auto work = [task](){ (*task)(); };
    unsigned int i = index_++;
    for(unsigned int n = 0; n < count_; ++n) {
      if(queues_[(i + n) % count_].TryPush(work)){
        return res;
      }
    }
    queues_[i % count_].Push(work);

    return res;
  }

  template<typename F, typename... Args>
  std::future<typename std::result_of<F(Args...)>::type> AsyncTask(F&& f, Args&&... args) {
    using return_type = typename std::result_of<F(Args...)>::type;
    auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    std::future<return_type> res = task->get_future();

    auto work = [task](){ (*task)(); };
    unsigned int i = index_++;
    for(unsigned int n = 0; n < count_; ++n) {
      if(queues_[(i + n) % count_].TryPush(work)){
        return res;
      }
    }
    queues_[i % count_].Push(work);

    return res;
  }

  unsigned int get_num_threads(){
    return count_;
  }

private:
  void Run(unsigned int i) {
    while (true) {
      std::function<void()> f;
      for (unsigned n = 0; n != count_; ++n) {
        if (queues_[(i + n) % count_].TryPop(f)) {
          break;
        }
      }
      if (!f && !queues_[i].Pop(f)) {
        break;
      }
      f();
    }
  }

  const unsigned int count_;
  std::vector<std::thread> threads_;
  std::atomic<unsigned int> index_;
  std::vector<NotificationQueue> queues_;
};

}
#endif
