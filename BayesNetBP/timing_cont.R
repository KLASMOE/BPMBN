#连续点比例变化

library(bnlearn)
library(BayesNetBP)
library(readxl)
library(igraph)
library(ggplot2)

#setwd("results")  #设置工作路径

index1 <- 1
i <- 1
tree.init.p <- list()
tree.post <- list()
marg <- list()

V <- 50 # 节点个数向量
C <- seq(0, 100, 25) # 连续点的比例（扩大了100倍）
E <- 5 # 证据变量个数
N <- 1000 # 样本量
P <- 0:99 # 第p个数据


time <- vector("numeric", length = length(C))

for (v in V) {
  file_name <- paste("Time_V", v, "E", E, ".txt", sep="") # 保存运行时间的文件
  if (file.exists(file_name)) {
    file.remove(file_name)  # 如果文件存在则删除
  }
  #cat(sprintf("V%iE%i\n", as.integer(v), as.integer(E)), file = file_name, append = TRUE)
  for (c in C) {
    # 定义文件名并初始化文件
    start_time <- Sys.time()
    for (p in P) {
      # 提取结构数据
      lines1 <- readLines(paste("../BPMBN/BPMBN/results/structure/structure_V", v, "C", c, "_", p, ".txt", sep=""))
      # 提取第2行到第51行的数据
      selected_lines <- lines1[2:(v+1)]

      # 找到第2行到第v行的最大列数
      max_cols <- max(sapply(selected_lines, function(line1) length(strsplit(line1, " ")[[1]])))

      # 导入数据，按照最大列数导入
      variable <- read.table(text = selected_lines, fill = TRUE, col.names = paste0("V", 1:max_cols))
      variable <- variable[,-2]
      variable1 <- variable[,1]

      # 将0和NA值不做处理
      variable[] <- lapply(variable, function(x) {
        ifelse(!is.na(x) & x != 0, paste0("X", x), x)
      })

      variable <- as.data.frame(lapply(variable[,2:(max_cols-1)], as.character), stringsAsFactors = FALSE)

      nodes <- paste0("X",1:v)
      net <- empty.graph(nodes)

      # 针对数据框的每一列，应用自定义函数
      for (j1 in 2:ncol(variable)) {
        for (i1 in 1:nrow(variable)){
          if (is.na(variable[i1,j1]) || variable[i1,j1] == 0) {
            next  # 如果值为NA或0，则跳转到下一行
          } else {
            # 进行某种运算
            net <- set.arc(net, variable[i1, j1], variable[i1, 1])
          }
        }
      }

      amat <- amat(net)
      # 转换为igraph图对象
      h <- graph_from_adjacency_matrix(amat, mode = "directed", diag = FALSE)

      # 设置输出路径
      output_path <- "./"

      # pdf(paste(output_path, "plot", v, ".pdf", sep=""))

      plot(h, layout=layout_nicely(h),
           vertex.size=6,
           vertex.label.color="black",
           vertex.color="skyblue",
           edge.arrow.size=0.2,
           edge.color = "red",  # 设置箭头颜色为红色
           main=paste("Bayes Network Structure V=", v ),
           vertex.label.cex=0.4)
      # 关闭PDF文件
      dev.off()

      # 计算向量中-1和-2的个数
      n1 <- sum(variable1 == -1)
      n2 <- sum(variable1 == -2)

      # 定义节点是离散节点（TRUE）还是连续节点（FALSE）
      node.class <- rep(c(TRUE, FALSE), c(n1,n2))
      names(node.class) <- paste0("X",1:v)
      node.class

      # 转换网络格式
      g <- as.graphNEL(net)

      for (e in E) {
        if (e==0) {
          for (n in N) {
            # 读取数据
            data_path <- paste0("../BPMBN/BPMBN/results/data/data_V", v, "C", c, "E", e, "N", n,  "_", p, ".csv")

            # 读取Excel文件
            data1 <- read.csv(data_path)
            # 将矩阵转换为data.frame
            # 为data.frame的每列变量命名为X1, X2, ..., Xk

            data_df <- as.data.frame(data1)
            names(data_df) <- paste0("X",1:v)

            # 尝试执行计算
            result1 <- tryCatch({
              # 计算结果
              Initializer(dag=g, data=data_df, node.class=node.class)
            }, error = function(e1) {
              # 如果出现错误，将结果记为NA
              NA
            })

            tree.init.p[[index1]] <- result1


            # 进行传播计算
            result3 <- tryCatch({
              # 进行传播计算，计算边际分布
              Marginals(tree.init.p[[index1]], names(data_df))
            }, error = function(e3) {
              # 如果出现错误，将结果记为NA
              NA
            })
            marg[[index1]] <- result3

            index1 <- index1 + 1
          }
        } else {
          # 读取证据
          ev <- read.csv(paste0("../BPMBN/BPMBN/results/evidence/evidence_V", v, "C", c, "E", e, ".csv"))
          # 获取列名
          ev1 <- colnames(ev)
          # 获取第一行数据
          ev2 <- unname(as.list(ev[1, ]))

          not_in_ev1 <- nodes[!(nodes %in% ev1)] # 查找所有非证据变量。

          for (n in N) {
            # 读取数据
            data_path <- paste0("../BPMBN/BPMBN/results/data/data_V", v, "C", c, "E", e, "N", n,  "_", p, ".csv")

            # 读取Excel文件
            data1 <- read.csv(data_path)
            # 将矩阵转换为data.frame
            # 为data.frame的每列变量命名为X1, X2, ..., Xk

            data_df <- as.data.frame(data1)
            names(data_df) <- paste0("X",1:v)

            # 尝试执行计算
            result1 <- tryCatch({
              # 计算结果
              Initializer(dag=g, data=data_df, node.class=node.class)
            }, error = function(e1) {
              # 如果出现错误，将结果记为NA
              NA
            })

            tree.init.p[[index1]] <- result1

            # 指定证据变量及其取值
            result2 <- tryCatch({
              # 设定证据变量及其值
              AbsorbEvidence(tree.init.p[[index1]], vars = ev1, values = ev2 )
            }, error = function(e2) {
              # 如果出现错误，将结果记为NA
              NA
            })
            tree.post[[index1]] <- result2


            # 进行传播计算
            result3 <- tryCatch({
              # 进行传播计算，计算边际分布
              Marginals(tree.post[[index1]], not_in_ev1)
            }, error = function(e3) {
              # 如果出现错误，将结果记为NA
              NA
            })
            marg[[index1]] <- result3

            index1 <- index1 + 1
          }
        }
      }
    }
    end_time <- Sys.time()
    time[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
    Average_Time <- time[i]/length(P)
    # cat(sprintf("V=%i, N=%i, E=%i, C=%.2f,: Times = %s (s), Average Time = %s (s) \n",
    #            as.integer(v), as.integer(n), as.integer(e), c/100,
    #            paste(round(time[i],5), collapse = ", "), round(Average_Time,5)), file = file_name, append = TRUE)
    # 假设在适当的循环或逻辑中
    cat(sprintf("%s ", round(Average_Time,5)), file = file_name, append = TRUE)

    i <- i + 1
  }
  mean_time <- time/length(P)
}





data <- data.frame(x = (C/100), y = mean_time)

# 使用 ggplot 绘制时间折线图，添加颜色映射以生成图例
ggplot(data, aes(x = x, y = y)) +
  geom_line(color = "blue") +  # 绘制线条
  geom_point(shape = 21, size = 3, fill = "red", color = "black") +  # 更改点的样式 +  # 添加点
  geom_text(aes(label = sprintf("%.3f", y)), nudge_y = 0.05) +  # 添加文本标签，使用 sprintf 格式化为三位小数
  labs(title = paste0("传播计算平均时间随连续点比例变化折线（V=", V, ",E=", E, ",N=", N, "）"), x = "连续点比例", y = "平均时间(seconds)") +  # 添加标签和图例标题
  theme_minimal() +  # 使用简洁主题
  theme(plot.title = element_text(hjust = 0.5))  # 确保标题居中（通常默认居中）

