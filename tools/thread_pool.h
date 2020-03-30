//
//  thread_pool.h
//  ppCLASS
//
//  Created by Thomas Tram on 02/03/2020.
//  Copyright © 2020 Aarhus University. All rights reserved.
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
  TaskSystem(unsigned int count = std::thread::hardware_concurrency())
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
  std::vector<NotificationQueue> queues_;
  std::atomic<unsigned int> index_;
};

}
#endif //THREAD_POOL_H
