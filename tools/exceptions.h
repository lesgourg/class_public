//
//  exceptions.h
//  ppCLASS
//
//  Created by Thomas Tram on 16/03/2020.
//  Copyright © 2020 Aarhus University. All rights reserved.
//

#ifndef exceptions_h
#define exceptions_h

#include <string>

void ThrowInvalidArgumentIf(bool condition, std::string string_for_printf, ...);

void ThrowInvalidArgument(std::string string_for_printf, ...);

#endif /* exceptions_h */
