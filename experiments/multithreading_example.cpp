#include <iostream>
#include <thread>
#include <vector>

void function()
{
    std::cout << "Tread ID = " << std::this_thread::get_id() << std::endl;
}

int main()
{
    int num_threads = 500; // Change this to the desired number of threads

    std::vector<std::thread> threads;

    // Create multiple threads
    for (int i = 0; i < num_threads; i++)
    {
        threads.push_back(std::thread(function));
    }

    // Wait for the threads to finish
    for (int i = 0; i < num_threads; i++)
    {
        threads[i].join();
    }

    std::cout << "All threads have finished!" << std::endl;

    return 0;
}