#!MC 1410
$!OpenLayout  "/home/geve/Dropbox/EnKF_scalar/Run/RunReview/review.lay"

$!ActiveLineMaps -= [1-12]
$!ActiveLineMaps += [1-4,11]
$!Redraw 
$!PrintSetup Palette = Color
$!ExportSetup ExportFormat = EPS
$!ExportSetup ImageWidth = 1569
$!ExportSetup EPSPreviewImage{ImageType = None}
$!ExportSetup ExportFName = '/home/geve/Dropbox/EnKF_scalar/Run/RunReview/scalarES.eps'
$!Export 
  ExportRegion = AllFrames

$!ActiveLineMaps -= [1-12]
$!ActiveLineMaps += [1-2,5-6]
$!Redraw 
$!PrintSetup Palette = Color
$!ExportSetup ExportFormat = EPS
$!ExportSetup ImageWidth = 1569
$!ExportSetup EPSPreviewImage{ImageType = None}
$!ExportSetup ExportFName = '/home/geve/Dropbox/EnKF_scalar/Run/RunReview/scalarIES.eps'
$!Export 
  ExportRegion = AllFrames

$!ActiveLineMaps -= [1-12]
$!ActiveLineMaps += [1-2,7-8,12]
$!Redraw 
$!PrintSetup Palette = Color
$!ExportSetup ExportFormat = EPS
$!ExportSetup ImageWidth = 1569
$!ExportSetup EPSPreviewImage{ImageType = None}
$!ExportSetup ExportFName = '/home/geve/Dropbox/EnKF_scalar/Run/RunReview/scalarESMDA.eps'
$!Export 
  ExportRegion = AllFrames

$!ActiveLineMaps -= [1-12]
$!ActiveLineMaps += [1-2,9-10]
$!Redraw 
$!PrintSetup Palette = Color
$!ExportSetup ExportFormat = EPS
$!ExportSetup ImageWidth = 1569
$!ExportSetup EPSPreviewImage{ImageType = None}
$!ExportSetup ExportFName = '/home/geve/Dropbox/EnKF_scalar/Run/RunReview/scalarSTEIN.eps'
$!Export 
  ExportRegion = AllFrames
