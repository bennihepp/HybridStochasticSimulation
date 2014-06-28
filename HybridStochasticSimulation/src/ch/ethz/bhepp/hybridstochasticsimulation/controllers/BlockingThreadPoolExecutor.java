package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

public class BlockingThreadPoolExecutor implements ExecutorService {

	public static ExecutorService blockingDecorator(ExecutorService executor, int bound) {
		return new BlockingThreadPoolExecutor(executor, bound);
	}

	private final ExecutorService executor;
	private final Semaphore semaphore;

	private BlockingThreadPoolExecutor(ExecutorService executor, int bound) {
		this.executor = executor;
		this.semaphore = new Semaphore(bound);
	}

	@Override
	public void execute(final Runnable command) {
		try {
			semaphore.acquire();
			try {
	            executor.execute(new Runnable() {

					@Override
					public void run() {
						try {
							command.run();
						} finally {
							semaphore.release();
						}
					}

	            });
	        } catch (Exception e) {
	        	// This also handles RejectedExecutionException
	        	semaphore.release();
	        	throw e;
	        }
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new RejectedExecutionException("Interrupted while trying to execute task", e);
		}
	}

	@Override
	public void shutdown() {
		executor.shutdown();
	}

	@Override
	public List<Runnable> shutdownNow() {
		return executor.shutdownNow();
	}

	@Override
	public boolean isShutdown() {
		return executor.isShutdown();
	}

	@Override
	public boolean isTerminated() {
		return executor.isTerminated();
	}

	@Override
	public boolean awaitTermination(long timeout, TimeUnit unit) throws InterruptedException {
		return executor.awaitTermination(timeout, unit);
	}

	@Override
	public <T> Future<T> submit(Callable<T> task) {
		return executor.submit(task);
	}

	@Override
	public <T> Future<T> submit(Runnable task, T result) {
		return executor.submit(task, result);
	}

	@Override
	public Future<?> submit(Runnable task) {
		return executor.submit(task);
	}

	@Override
	public <T> List<Future<T>> invokeAll(Collection<? extends Callable<T>> tasks) throws InterruptedException {
		return executor.invokeAll(tasks);
	}

	@Override
	public <T> List<Future<T>> invokeAll(Collection<? extends Callable<T>> tasks, long timeout, TimeUnit unit)
			throws InterruptedException {
		return executor.invokeAll(tasks, timeout, unit);
	}

	@Override
	public <T> T invokeAny(Collection<? extends Callable<T>> tasks) throws InterruptedException, ExecutionException {
		return executor.invokeAny(tasks);
	}

	@Override
	public <T> T invokeAny(Collection<? extends Callable<T>> tasks, long timeout, TimeUnit unit)
			throws InterruptedException, ExecutionException, TimeoutException {
		return executor.invokeAny(tasks, timeout, unit);
	}

}
