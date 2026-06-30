path1 = '.\\mPFCm9\\'
path2 = '_Traces.csv'

BLpath = paste(path1,'BL',path2, sep = '')
BLraw = read.csv(BLpath)
BL = BLraw[-c(1),-c(1)]
BL = sapply(BL,as.numeric)

SDpath = paste(path1,'SD',path2, sep = '')
SDraw = read.csv(SDpath)
SD = SDraw[-c(1),-c(1)]

WOpath = paste(path1,'WO',path2, sep = '')
WOraw = read.csv(WOpath)
WO = WOraw[-c(1),-c(1)]
