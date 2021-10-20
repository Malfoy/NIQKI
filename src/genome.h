#ifndef __GENOME_H__
#define __GENOME_H__

#include <string>


  class Genome {

  public:


    inline void setName(const std::string &genome_name) {
      name = genome_name;
    }


    inline std::string getName() const {
      return name;
    }

    /**
     * \brief Set the index filename.
     *
     * \param idx The index filename to set.
     */
    inline void setIndexFile(const std::string &idx) {
      index_file = idx;
    }

    /**
     * \brief Returns the index filename if any.
     *
     * \return Returns the index filename if any.
     */
    inline std::string getIndexFile() const {
      return index_file;
    }

  private:

    std::string name;  //genome name
    std::string index_file; // if load from index file

};

#endif
