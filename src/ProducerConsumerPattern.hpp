/*
 * ProducerConsumerPattern.h
 *
 *  Created on: Feb 1, 2017
 *      Author: sarah
 */

#pragma once

#include <stddef.h>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <vector>
#include <iostream>
#include <cassert>
#include <thread>

template<class T> class ProducerConsumerPattern {
public:
	ProducerConsumerPattern(size_t maxDataQueueSize, std::function<double(std::vector<T>&, size_t)> produce, std::function<void(std::vector<T>&, size_t)> consume);
	void run(size_t numProducers, size_t numConsumers);
private:
	void produce(size_t producerId);
	void consume(size_t consumerId, size_t myProducer);

	std::function<double(std::vector<T>&, size_t)> prepareData;
	std::function<void(std::vector<T>&, size_t)> processData;

	std::vector<std::mutex> mtx;
	std::vector<std::condition_variable> cvConsume;
	std::vector<std::condition_variable> cvProduce;
	std::vector<std::queue<std::vector<T> > > dataQueue;

	std::vector<bool> done;

	size_t _maxDataQueueSize;
};


template<class T> ProducerConsumerPattern<T>::ProducerConsumerPattern(size_t maxDataQueueSize,
		std::function<double(std::vector<T>&, size_t)> produce,
		std::function<void(std::vector<T>&, size_t)> consume) {
	_maxDataQueueSize = maxDataQueueSize;
	prepareData = produce;
	processData = consume;
}

// assumes that
template<class T> void ProducerConsumerPattern<T>::run(size_t numProducers, size_t numConsumers) {
	std::vector<std::thread> threads;

	cvConsume = std::vector<std::condition_variable>(numProducers);
	cvProduce = std::vector<std::condition_variable>(numProducers);
	mtx = std::vector<std::mutex>(numProducers);
	dataQueue = std::vector<std::queue<std::vector<T> > >(numProducers);
	done = std::vector<bool>(numProducers, false);

	assert(numConsumers % numProducers == 0);
	size_t consumersPerProducer = numConsumers / numProducers;

	for (size_t i = 0; i < numProducers; ++i) {
		threads.push_back(std::thread(&ProducerConsumerPattern<T>::produce, this, i));
	}
	for (size_t i = 0; i < numConsumers; ++i) {
		threads.push_back(std::thread(&ProducerConsumerPattern<T>::consume, this, i, i / consumersPerProducer));
	}
	for (size_t i = 0; i < threads.size(); ++i) {
		threads[i].join();
	}
	std::cout << "Threads finished.\n";
}

template<class T> void ProducerConsumerPattern<T>::produce(size_t producerId) {
	size_t actReads = 0;
	double minProgress = 1;
	while (!done[producerId]) {
		std::vector<T> buffer;
		double progress = prepareData(buffer, producerId);
		actReads += buffer.size();
		std::unique_lock<std::mutex> lck(mtx[producerId]);

		cvProduce[producerId].wait(lck, [this, &producerId] {return dataQueue[producerId].size() <= _maxDataQueueSize;});
		dataQueue[producerId].push(buffer);
		cvConsume[producerId].notify_one();

		if (progress >= minProgress) {
			std::cout << "Producer " + std::to_string(producerId) + ": " + std::to_string(progress) + " \%\n";
			minProgress += 1;
		}
		if (progress == 100) {
			done[producerId] = true;
			cvConsume[producerId].notify_all();
		}
	}
	std::cout << "Producer " + std::to_string(producerId) + ": finished producing.\n";
}

template<class T> void ProducerConsumerPattern<T>::consume(size_t consumerId, size_t myProducer) {
	while (!done[myProducer] || !dataQueue.empty()) {
		std::unique_lock<std::mutex> lck(mtx[myProducer]);
		cvConsume[myProducer].wait(lck, [this, &myProducer] {return (!dataQueue[myProducer].empty()) || done[myProducer];});
		{
			if (dataQueue[myProducer].empty()) {
				break;
			}
			std::vector<T> buffer = dataQueue[myProducer].front();
			dataQueue[myProducer].pop();
			lck.unlock();
			cvProduce[myProducer].notify_one();
			processData(buffer, consumerId);
		}
	}
	std::cout << "Consumer " + std::to_string(consumerId) + ": finished consuming.\n";
}
