/*!
 * @file thread_pool.hpp
 *
 * @brief ThreadPool class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>
#include <queue>
#include <mutex>
#include <thread>
#include <future>
#include <atomic>
#include <condition_variable>

class Semaphore;
class ThreadPool;

std::unique_ptr<Semaphore> createSemaphore(uint32_t value);

std::unique_ptr<ThreadPool> createThreadPool(size_t num_threads =
    std::thread::hardware_concurrency());

class Semaphore {
public:

    ~Semaphore() = default;

    uint32_t value() const {
        return value_;
    }

    void wait();
    void post();

    friend std::unique_ptr<Semaphore> createSemaphore(uint32_t value);

private:

    Semaphore(uint32_t value);
    Semaphore(const Semaphore&) = delete;
    const Semaphore& operator=(const Semaphore&) = delete;

    std::mutex mutex_;
    std::condition_variable condition_;
    uint32_t value_;
};

class ThreadPool {
public:

    ~ThreadPool();

    size_t num_threads() const {
        return threads_.size();
    }

    template<typename T, typename... Ts>
    auto submit_task(T&& routine, Ts&&... params)
        -> std::future<typename std::result_of<T(Ts...)>::type> {

        auto task = std::make_shared<std::packaged_task<typename std::result_of<T(Ts...)>::type()>>(
            std::bind(std::forward<T>(routine), std::forward<Ts>(params)...)
        );

        auto task_result = task->get_future();

        auto task_wrapper = [task]() {
            (*task)();
        };

        queue_sem_->wait();

        task_queue_.emplace(task_wrapper);

        queue_sem_->post();
        active_sem_->post();

        return task_result;
    }

    friend std::unique_ptr<ThreadPool> createThreadPool(size_t num_threads);

private:

    ThreadPool(size_t num_threads);
    ThreadPool(const ThreadPool&) = delete;
    const ThreadPool& operator=(const ThreadPool&) = delete;

    static void worker_thread(ThreadPool* thread_pool);

    std::vector<std::thread> threads_;

    std::queue<std::function<void()>> task_queue_;

    std::unique_ptr<Semaphore> queue_sem_;
    std::unique_ptr<Semaphore> active_sem_;

    std::atomic<bool> terminate_;
};
