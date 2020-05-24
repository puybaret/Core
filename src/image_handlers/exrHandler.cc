/****************************************************************************
 *
 * 		exrHandler.cc: EXR format handler
 *      This is part of the yafray package
 *      Copyright (C) 2010 Rodrigo Placencia Vazquez
 *
 *      This library is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public
 *      License as published by the Free Software Foundation; either
 *      version 2.1 of the License, or (at your option) any later version.
 *
 *      This library is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *      Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this library; if not, write to the Free Software
 *      Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include <core_api/environment.h>
#include <core_api/imagehandler.h>
#include <core_api/params.h>
#include <core_api/scene.h>
#include <utilities/math_utils.h>
#include <utilities/fileUtils.h>

#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfRgbaFile.h>
#include <ImfArray.h>
#include <ImfVersion.h>

using namespace Imf;
using namespace Imath;

__BEGIN_YAFRAY

typedef genericScanlineBuffer_t<Rgba> halfRgbaScanlineImage_t;
typedef genericScanlineBuffer_t<float> grayScanlineImage_t;

class exrHandler_t: public imageHandler_t
{
public:
	exrHandler_t();
	~exrHandler_t();
	bool loadFromFile(const std::string &name);
	bool saveToFile(const std::string &name, int imgIndex = 0);
    bool saveToFileMultiChannel(const std::string &name, const renderPasses_t *renderPasses);
	static imageHandler_t *factory(paraMap_t &params, renderEnvironment_t &render);
	bool isHDR() { return true; }
};

exrHandler_t::exrHandler_t()
{
	handlerName = "EXRHandler";
}

exrHandler_t::~exrHandler_t()
{
	clearImgBuffers();
}

bool exrHandler_t::saveToFile(const std::string &name, int imgIndex)
{
	int h = getHeight(imgIndex);
	int w = getWidth(imgIndex);

	std::string nameWithoutTmp = name;
	nameWithoutTmp.erase(nameWithoutTmp.length()-4);
	if(session.renderInProgress()) Y_INFO << handlerName << ": Autosaving partial render (" << RoundFloatPrecision(session.currentPassPercent(), 0.01) << "% of pass " << session.currentPass() << " of " << session.totalPasses() << ") RGB" << ( m_hasAlpha ? "A" : "" ) << " file as \"" << nameWithoutTmp << "\"...  " << getDenoiseParams()  << yendl;
	else Y_INFO << handlerName << ": Saving RGB" << ( m_hasAlpha ? "A" : "" ) << " file as \"" << nameWithoutTmp << "\"...  " << getDenoiseParams()  << yendl;

	int chan_size = sizeof(half);
	const int num_colchan = 4;
	int totchan_size = num_colchan*chan_size;

	Header header(w, h);

	header.compression() = ZIP_COMPRESSION;

	header.channels().insert("R", Channel(HALF));
	header.channels().insert("G", Channel(HALF));
	header.channels().insert("B", Channel(HALF));
	header.channels().insert("A", Channel(HALF));

	OutputFile file(name.c_str(), header);

	Imf::Array2D<Imf::Rgba> pixels;
	pixels.resizeErase(h, w);

	for(int i = 0; i < w; ++i)
	{
		for(int j = 0; j < h; ++j)
		{
			colorA_t col = imgBuffer.at(imgIndex)->getColor(i, j);
			pixels[j][i].r = col.R;
			pixels[j][i].g = col.G;
			pixels[j][i].b = col.B;
			pixels[j][i].a = col.A;
		}
	}

	char* data_ptr = (char *)&pixels[0][0];
	
	FrameBuffer fb;
	fb.insert("R", Slice(HALF, data_ptr              , totchan_size, w * totchan_size));
	fb.insert("G", Slice(HALF, data_ptr +   chan_size, totchan_size, w * totchan_size));
	fb.insert("B", Slice(HALF, data_ptr + 2*chan_size, totchan_size, w * totchan_size));
	fb.insert("A", Slice(HALF, data_ptr + 3*chan_size, totchan_size, w * totchan_size));

	file.setFrameBuffer(fb);

	try
	{
		file.writePixels(h);
		Y_VERBOSE << handlerName << ": Done." << yendl;
		return true;
	}
	catch (const std::exception &exc)
	{
		Y_ERROR << handlerName << ": " << exc.what() << yendl;
		return false;
	}
}

bool exrHandler_t::saveToFileMultiChannel(const std::string &name, const renderPasses_t *renderPasses)
{
	int h0 = imgBuffer.at(0)->getHeight();
	int w0 = imgBuffer.at(0)->getWidth();

	bool allImageBuffersSameSize = true;
	for(size_t idx = 0; idx < imgBuffer.size(); ++idx)
	{
		if(imgBuffer.at(idx)->getHeight() != h0) allImageBuffersSameSize = false;
		if(imgBuffer.at(idx)->getWidth() != w0) allImageBuffersSameSize = false;
	}
	
	if(!allImageBuffersSameSize)
	{
		Y_ERROR << handlerName << ": Saving Multilayer EXR failed: not all the images in the imageBuffer have the same size. Make sure all images in buffer have the same size or use a non-multilayered EXR format." << yendl;
		return false;
	}

    std::string extPassName;

	std::string nameWithoutTmp = name;
	nameWithoutTmp.erase(nameWithoutTmp.length()-4);
    
    if(session.renderInProgress()) Y_INFO << handlerName << ": Autosaving partial render (" << RoundFloatPrecision(session.currentPassPercent(), 0.01) << "% of pass " << session.currentPass() << " of " << session.totalPasses() << ") Multilayer EXR" << " file as \"" << nameWithoutTmp << "\"...  " << getDenoiseParams() << yendl;
    else Y_INFO << handlerName << ": Saving Multilayer EXR" << " file as \"" << nameWithoutTmp << "\"...  " << getDenoiseParams()  << yendl;

	int chan_size = sizeof(half);
	const int num_colchan = 4;
	int totchan_size = num_colchan*chan_size;
	
    Header header(w0, h0);
    FrameBuffer fb;
	header.compression() = ZIP_COMPRESSION;
    
	std::vector<Imf::Array2D<Imf::Rgba> *> pixels;

    for(size_t idx = 0; idx < imgBuffer.size(); ++idx)
    {
		extPassName = "RenderLayer." + renderPasses->extPassTypeStringFromIndex(idx) + ".";        
		Y_VERBOSE << "    Writing EXR Layer: " << renderPasses->extPassTypeStringFromIndex(idx) << yendl;
        
        const std::string channelR_string = extPassName + "R";
        const std::string channelG_string = extPassName + "G";
        const std::string channelB_string = extPassName + "B";
        const std::string channelA_string = extPassName + "A";
        
        const char* channelR = channelR_string.c_str();
        const char* channelG = channelG_string.c_str();
        const char* channelB = channelB_string.c_str();
        const char* channelA = channelA_string.c_str();
        
        header.channels().insert(channelR, Channel(HALF));
        header.channels().insert(channelG, Channel(HALF));
        header.channels().insert(channelB, Channel(HALF));
        header.channels().insert(channelA, Channel(HALF));
 
		pixels.push_back(new Imf::Array2D<Imf::Rgba>);
		pixels.at(idx)->resizeErase(h0, w0);

		for(int i = 0; i < w0; ++i)
		{
			for(int j = 0; j < h0; ++j)
			{
				colorA_t col = imgBuffer.at(idx)->getColor(i, j);
				(*pixels.at(idx))[j][i].r = col.R;
				(*pixels.at(idx))[j][i].g = col.G;
				(*pixels.at(idx))[j][i].b = col.B;
				(*pixels.at(idx))[j][i].a = col.A;
			}
		}

		char* data_ptr = (char *)&(*pixels.at(idx))[0][0];

        fb.insert(channelR, Slice(HALF, data_ptr              , totchan_size, w0 * totchan_size));
        fb.insert(channelG, Slice(HALF, data_ptr +   chan_size, totchan_size, w0 * totchan_size));
        fb.insert(channelB, Slice(HALF, data_ptr + 2*chan_size, totchan_size, w0 * totchan_size));
        fb.insert(channelA, Slice(HALF, data_ptr + 3*chan_size, totchan_size, w0 * totchan_size));
    }
    
    OutputFile file(name.c_str(), header);    
	file.setFrameBuffer(fb);
	
	try
	{
		file.writePixels(h0);
		Y_VERBOSE << handlerName << ": Done." << yendl;
		for(size_t idx = 0; idx < pixels.size(); ++idx)
		{
			delete pixels.at(idx);
			pixels.at(idx) = nullptr;
		}
		pixels.clear();
		return true;
	}
	catch (const std::exception &exc)
	{
		Y_ERROR << handlerName << ": " << exc.what() << yendl;
		for(size_t idx = 0; idx < pixels.size(); ++idx)
		{
			delete pixels.at(idx);
			pixels.at(idx) = nullptr;
		}
		pixels.clear();
		return false;
	}
}

bool exrHandler_t::loadFromFile(const std::string &name)
{
	std::string tempFilePathString = "";    //filename of the temporary exr file that will be generated to deal with the UTF16 ifstream path problems in OpenEXR libraries with MinGW

	FILE *fp = fileUnicodeOpen(name.c_str(), "rb");
	Y_INFO << handlerName << ": Loading image \"" << name << "\"..." << yendl;
	
	if(!fp)
	{
		Y_ERROR << handlerName << ": Cannot open file " << name << yendl;
		return false;
	}
	else
	{
		char bytes[4];
		fread(&bytes, 1, 4, fp);
#if defined(_WIN32)		
		fseek (fp , 0 , SEEK_SET);
		auto tempFolder = boost::filesystem::temp_directory_path();
		auto tempFile = boost::filesystem::unique_path();
		tempFilePathString = tempFolder.string() + tempFile.string() + ".exr";
		Y_VERBOSE << handlerName << ": Creating intermediate temporary file tempstr=" << tempFilePathString << yendl;
		FILE *fpTemp = fopen(tempFilePathString.c_str(), "wb");
		if(!fpTemp)
		{
			Y_ERROR << handlerName << ": Cannot create intermediate temporary file " << tempFilePathString << yendl;
			return false;
		}
		//Copy original EXR texture contents into new temporary file, so we can circumvent the lack of UTF16 support in MinGW ifstream
		unsigned char *copy_buffer = new unsigned char [1024];
		int numReadBytes = 0;
		while((numReadBytes = fread(copy_buffer, sizeof(unsigned char), 1024, fp)) == 1024) fwrite(copy_buffer, sizeof(unsigned char), 1024, fpTemp);
		fwrite(copy_buffer, sizeof(unsigned char), numReadBytes, fpTemp);
		fclose(fpTemp);		
		delete [] copy_buffer;
#endif		
		fclose(fp);
		fp = nullptr;
		if(!isImfMagic(bytes)) return false;
	}

	try
	{
#if defined(_WIN32)
		Y_INFO << handlerName << ": Loading intermediate temporary file tempstr=" << tempFilePathString << yendl;
		RgbaInputFile file(tempFilePathString.c_str());		
#else		
		RgbaInputFile file(name.c_str());		
#endif		
		Box2i dw = file.dataWindow();

		m_width  = dw.max.x - dw.min.x + 1;
		m_height = dw.max.y - dw.min.y + 1;
		m_hasAlpha = true;

		clearImgBuffers();

		int nChannels = 3;
		if(m_grayscale) nChannels = 1;
		else if(m_hasAlpha) nChannels = 4;

		imgBuffer.push_back(new imageBuffer_t(m_width, m_height, nChannels, getTextureOptimization()));

		Imf::Array2D<Imf::Rgba> pixels;
		pixels.resizeErase(m_width, m_height);
		file.setFrameBuffer(&pixels[0][0] - dw.min.y - dw.min.x * m_height, m_height, 1);
		file.readPixels(dw.min.y, dw.max.y);
		
		for(int i = 0; i < m_width; ++i)
		{
			for(int j = 0; j < m_height; ++j)
			{
				colorA_t col;
				col.R = pixels[i][j].r;
				col.G = pixels[i][j].g;
				col.B = pixels[i][j].b;
				col.A = pixels[i][j].a;
				imgBuffer.at(0)->setColor(i, j, col, m_colorSpace, m_gamma);
			}
		}
	}
	catch (const std::exception &exc)
	{
		Y_ERROR << handlerName << ": " << exc.what() << yendl;
		return false;
	}

#if defined(_WIN32)
	Y_INFO << handlerName << ": Deleting intermediate temporary file tempstr=" << tempFilePathString << yendl;
	std::remove(tempFilePathString.c_str());		
#endif	

	return true;
}

imageHandler_t *exrHandler_t::factory(paraMap_t &params,renderEnvironment_t &render)
{
	int pixtype = HALF;
	int compression = ZIP_COMPRESSION;
	int width = 0;
	int height = 0;
	bool withAlpha = false;
	bool forOutput = true;
	bool multiLayer = false;
	bool img_grayscale = false;
	bool denoiseEnabled = false;
	int denoiseHLum = 3;
	int denoiseHCol = 3;
	float denoiseMix = 0.8f;

	params.getParam("pixel_type", pixtype);
	params.getParam("compression", compression);
	params.getParam("width", width);
	params.getParam("height", height);
	params.getParam("alpha_channel", withAlpha);
	params.getParam("for_output", forOutput);
	params.getParam("img_multilayer", multiLayer);
	params.getParam("img_grayscale", img_grayscale);
/*	//Denoise is not available for HDR/EXR images
 * 	params.getParam("denoiseEnabled", denoiseEnabled);
 *	params.getParam("denoiseHLum", denoiseHLum);
 *	params.getParam("denoiseHCol", denoiseHCol);
 *	params.getParam("denoiseMix", denoiseMix);
 */
	imageHandler_t *ih = new exrHandler_t();
	
	ih->setTextureOptimization(TEX_OPTIMIZATION_HALF_FLOAT);

	if(forOutput)
	{
		if(yafLog.getUseParamsBadge()) height += yafLog.getBadgeHeight();
		ih->initForOutput(width, height, render.getRenderPasses(), denoiseEnabled, denoiseHLum, denoiseHCol, denoiseMix, withAlpha, multiLayer, img_grayscale);
	}

	return ih;
}

extern "C"
{

	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerImageHandler("exr", "exr", "EXR [IL&M OpenEXR]", exrHandler_t::factory);
	}

}
__END_YAFRAY
