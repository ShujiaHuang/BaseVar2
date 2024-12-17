## Learn from Seqlib 

https://github.com/jarro2783/cxxopts

## About ThreadPool

除了本程序在使用的 <threadpool.h> 之外，[这篇博客](https://murphypei.github.io/blog/2019/04/cpp-concurrent-4) 里还提到 github 上的另一个 ThreadPool: <https://github.com/mtrebi/thread-pool> 看起来似乎更好，但我还没有测试过，先记下。

- [2024-08-17 11:35:22] 它这个文档写得有点意思，会注明每个函数要干什么，以及为什么要做这个事情，可以参考用来完善 ilus。 

## Taskflow 

这个代码库可以用来更好地实现任务的并行和多线程：
- <https://www.zywvvd.com/notes/coding/cpp/taskflow/taskflow/>
- <https://www.cnblogs.com/ljmiao/p/18135164>

之后研究一下，看看是否有必要添加到 BaseVar2 中，替换掉原来的 threadpool ？

> **使用的时候，从 github 下载 taskflow，然后删掉目录下除了 taskflow 之外的所有文件和目录，即可使用**。
