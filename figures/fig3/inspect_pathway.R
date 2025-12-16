suppressPackageStartupMessages({library(CellChat)})
data('CellChatDB.human')
print(unique(CellChatDB.human$interaction$annotation))
