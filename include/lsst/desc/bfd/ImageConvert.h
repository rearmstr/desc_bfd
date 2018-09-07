#include "Image.h"
#include "lsst/afw/image/Image.h"
#include <boost/make_shared.hpp>


namespace lsst { namespace desc {  namespace bfd {

template<typename PixelT=float>
class ImageConverter
{
public:
    ImageConverter(PTR(afw::image::Image<PixelT>) image, afw::geom::Box2I box):
        _image(image), _box(box) {}
    ImageConverter(PTR(afw::image::Image<PixelT>) image):
        _image(image), _box(image->getBBox(afw::image::PARENT)) {}

    PTR(img::Image<PixelT>) getBFDImage() const {
        Bounds<int >const bounds(_box.getMinX(), _box.getMaxX(),
                                _box.getMinY(), _box.getMaxY());

        typename afw::image::Image<PixelT>::Array array = _image->getArray();
        int const stride = array.template getStride<0>();


        // Pointer to rows
        img::ImageData<PixelT> *img = new img::ImageData<PixelT>(bounds, false);
        img->setNotOwner();
        PixelT **row = img->pointer();

        // Need to account for the fact that dptr should not point to beginning of array
        PixelT *dptr = array.getData() - _box.getMinX();

        for(int i = _box.getMinY(); i <= _box.getMaxY(); ++i) {
            row[i] = dptr;
            dptr +=  stride;
        }
        return std::make_shared<img::Image<PixelT> >(img, new img::Header());
    }

private:
    PTR(afw::image::Image<PixelT>) _image;
    afw::geom::Box2I _box;
};

}}} // namespace lsst::meas::extensions::bfd
